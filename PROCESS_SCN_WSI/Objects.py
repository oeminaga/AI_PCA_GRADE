# %%
import xml.etree.ElementTree as ET
import slideio
import numpy as np
import os
import cv2
from skimage.measure import label, regionprops
from skimage.color import label2rgb
import random
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import math
import tensorflow as tf
from PIL import ImageFilter
import PIL

from PIL import Image


def GetMask(img, DebugMode=False):
    img_x = img.filter(ImageFilter.GaussianBlur(8))
    rgb = np.array(img_x.convert("RGB")).astype(np.uint8)
    hls = cv2.cvtColor(rgb, cv2.COLOR_RGB2HLS)
    mask = 255 - \
        cv2.inRange(hls, np.array([0, 225, 0]), np.array([255, 255, 255]))
    rgb = cv2.bitwise_and(rgb, rgb, mask=mask)
    # BLUE
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    lwr = np.array([90, 50, 50])
    upr = np.array([130, 255, 255])
    msk = cv2.inRange(hsv, lwr, upr)

    mask = mask-msk
    mask[mask < 0] = 0
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    hsv = cv2.bitwise_and(hsv, hsv, mask=mask)
    msk = cv2.inRange(hsv, (36, 0, 0), (90, 255, 255))

    mask = mask-msk
    mask[mask < 0] = 0

    kernel = np.ones((5, 5), np.uint8)
    mask = cv2.erode(mask, kernel, iterations=1)
    mask = cv2.dilate(mask, kernel, iterations=1)
    if DebugMode:
        plt.imshow(mask)
        plt.show()
    return mask


class HistoImage(tf.keras.utils.Sequence):
    def __init__(self, filename, batch_size=16, target_size=(512, 512), target_magnification=20, overlap_rate=(0.99,  0.99), preprocessing_function=None, verbose=0) -> None:
        super().__init__()
        self.filename = filename
        type_driver = filename.split(".")[-1]
        slide = slideio.open_slide(filename, type_driver.upper())
        self.image = slide.get_scene(1)
        self.filename_mask = self.filename.replace(f".{type_driver}", ".xml")
        self.ImageMagnification = int(self.image.magnification)
        self.batch_size = batch_size
        self.target_size = target_size
        self.overlap_rate = overlap_rate
        self.target_magnification = target_magnification
        self.preprocessing_function = preprocessing_function
        self.DebugMode = False
        # print(self.image.resolution)
        self.DebugModeV2 = False
        if verbose > 100:
            self.DebugMode = True

        if verbose > 1000:
            self.DebugModeV2 = True
        if os.path.exists(self.filename_mask):
            self.LoadPolygon(self.filename_mask)
            self.GenerateMask()
        self.GetThumbnails()

        self.GetImagePatchCoordinationFromMaskFile(
            target_magnification=self.target_magnification,
            target_patch_size=self.target_size,
            overlap_rate=self.overlap_rate,
            DebugMode=self.DebugMode,
            DebugModeV2=self.DebugModeV2
        )

    def extract_polygons(self, xml_root):
        annots = xml_root
        polys = []
        for annot in annots:
            for child in annot:
                if child.tag == 'Regions':
                    regions = child
                    for region in regions:
                        if region.tag == 'Region':
                            for child in region:
                                if child.tag == 'Vertices':
                                    vertices = child
                                    pts = []
                                    for vertex in vertices:
                                        pts.append(
                                            [int(vertex.attrib['X']), int(vertex.attrib['Y'])])

                                    pts = np.array([pts], np.int32)
                                    polys.append(pts)
        return polys

    def LoadPolygon(self, filename):
        self.polygons = self.extract_polygons(ET.parse(filename).getroot())

    def GetThumbnails(self):
        rect_1 = self.image.rect
        self.ThumbnailImage = Image.fromarray(self.image.read_block(
            (0, 0, rect_1[2], rect_1[3]), (rect_1[2]//100, rect_1[3]//100)))
        w, h = self.ThumbnailImage.size
        self.y_factor = 100  # self.image.dimensions[1]/h
        self.x_factor = 100  # self.image.dimensions[0]/w

    def GenerateMask(self):
        width, height = self.image.size
        dim = width, height
        self.mask = np.zeros(dim, dtype=np.uint8)
        self.dimension = (width, height)

        for polygon in self.polygons:
            self.mask = cv2.fillPoly(self.mask, polygon, (1))
        self.ThumbnailMask = Image.fromarray(cv2.resize(
            self.mask, (dim[1]//100, dim[0]//100))).rotate(90, expand=True)

    def GetImagePatchCoordinationFromMaskFile(self,
                                              MinFilledSize=0,  # 250,
                                              target_magnification=20,
                                              target_patch_size=(512, 512),
                                              overlap_rate=(0.99,  0.99),
                                              background_threshold=0.20,
                                              background_color=200,
                                              DebugMode=False,
                                              DebugModeV2=False,
                                              Limit=0):

        zoom = self.ImageMagnification
        w, h = self.dimension
        self.magnification_factor = target_magnification/zoom
        w_tiny, h_tiny = self.ThumbnailImage.size
        label_image = label(np.array(self.ThumbnailMask))
        regions_bb = []
        for region in regionprops(label_image):
            # take regions with large enough areas
            if region.filled_area >= MinFilledSize:
                regions_bb.append(region)
        if DebugMode:
            plt.imshow(self.ThumbnailMask)
            plt.imshow(self.ThumbnailImage, alpha=0.5)
            plt.show()
        store_coordination = []
        for region in regions_bb:
            reg = region.bbox
            min_row, min_col, max_row,  max_col = reg
            row_l_ = max_row - min_row
            col_l_ = max_col - min_col

            # x,y coordinate top_left
            Pt_top_left = ((min_col/w_tiny)*w*self.magnification_factor,
                           (min_row/h_tiny)*h*self.magnification_factor)
            # x,y coordinate bottom right
            Pt_bottom_right = ((max_col/w_tiny)*w*self.magnification_factor,
                               (max_row/h_tiny)*h*self.magnification_factor)
            org_x_l_ = Pt_bottom_right[0]-Pt_top_left[0]
            org_y_l_ = Pt_bottom_right[1]-Pt_top_left[1]
            number_of_patches_x_axis = org_x_l_/target_patch_size[0]
            number_of_patches_y_axis = org_y_l_/target_patch_size[1]

            patch_size_in_tiny_y_axis = round(row_l_/number_of_patches_y_axis)
            patch_size_in_tiny_x_axis = round(col_l_/number_of_patches_x_axis)

            if patch_size_in_tiny_y_axis != patch_size_in_tiny_x_axis:
                if patch_size_in_tiny_y_axis < patch_size_in_tiny_x_axis:
                    patch_size_in_tiny_x_axis = patch_size_in_tiny_y_axis
                else:
                    patch_size_in_tiny_y_axis = patch_size_in_tiny_x_axis

            masks = np.zeros((h_tiny, w_tiny))

            masks[min_row:max_row, min_col:max_col] = region.image

            thumbnail_img_corrected = np.array(self.ThumbnailImage).copy()
            for i in range(3):
                thumbnail_img_corrected[...,
                                        i] = thumbnail_img_corrected[..., i]*(masks)
            for y in range(min_row, max_row-math.floor(patch_size_in_tiny_y_axis*overlap_rate[1]), math.floor(patch_size_in_tiny_y_axis*overlap_rate[1])):
                for x in range(min_col, max_col-math.floor(patch_size_in_tiny_y_axis*overlap_rate[0]), math.floor(patch_size_in_tiny_x_axis*overlap_rate[0])):
                    ic = thumbnail_img_corrected[y:y +
                                                 patch_size_in_tiny_y_axis, x:x+patch_size_in_tiny_x_axis]
                    img_hsv = cv2.cvtColor(ic, cv2.COLOR_RGB2HSV)
                    bg_img = cv2.cvtColor(ic, cv2.COLOR_RGB2GRAY)
                    total_count = patch_size_in_tiny_x_axis*patch_size_in_tiny_y_axis
                    # lower mask (0-10)
                    lower_red = np.array([0, 50, 50])
                    upper_red = np.array([10, 255, 255])
                    mask0 = cv2.inRange(img_hsv, lower_red, upper_red)

                    # upper mask (170-180)
                    lower_red = np.array([170, 50, 50])
                    upper_red = np.array([180, 255, 255])
                    mask1 = cv2.inRange(img_hsv, lower_red, upper_red)
                    mask_b = mask0+mask1

                    mask_b[mask_b > 0] = 1
                    non_z_count_b = np.count_nonzero(mask_b)
                    per_red_background = non_z_count_b/total_count
                    if per_red_background > 0.5:
                        continue
                    
                    if DebugModeV2:
                        plt.imshow(bg_img, cmap="gray")
                        plt.show()
                        plt.imshow(bg_img>background_color, cmap="gray")
                        plt.show()
                    
                    non_z_count_ = np.count_nonzero(bg_img == 0)
                    non_z_count = np.count_nonzero(
                        bg_img > background_color)+non_z_count_
                    per_white_background = non_z_count/total_count

                    # if DebugMode and per_white_background<background_threshold:
                    #    ax.add_patch(Rectangle((x,y),patch_size_in_tiny_y_axis, patch_size_in_tiny_x_axis, edgecolor='r', facecolor="none"))
                    if per_white_background < background_threshold:
                        y_org_location = int(round(y*self.y_factor))
                        x_org_location = int(round(x*self.x_factor))
                        store_coordination.append([x_org_location, y_org_location, target_patch_size[0]*(
                            zoom//target_magnification), target_patch_size[1]*(zoom//target_magnification)])
        if DebugMode:
            fig, ax = plt.subplots()
            ax.imshow(self.ThumbnailImage)

            for coord in store_coordination:
                x = round(coord[0]/self.x_factor)
                y = round(coord[1]/self.y_factor)
                ax.add_patch(Rectangle((x, y), patch_size_in_tiny_y_axis,
                             patch_size_in_tiny_x_axis, edgecolor='r', facecolor="none"))
            plt.show()
        if Limit != 0:
            if len(store_coordination) > Limit:
                store_coordination = random.sample(store_coordination, Limit)
        self.store_coordination = store_coordination
        # return store_coordination

    def __len__(self):
        return round(len(self.store_coordination)/self.batch_size) if len(self.store_coordination) > self.batch_size else 1

    def __getitem__(self, idx):
        X_batch_tmp = self.store_coordination[idx *
                                              self.batch_size:(idx+1)*self.batch_size]
        X_batch = [np.array(Image.fromarray(self.image.read_block(
            x_b)).resize(self.target_size)) for x_b in X_batch_tmp]
        if self.preprocessing_function is not None:
            for i in range(X_batch.shape[0]):
                X_batch[i] = self.preprocessing_function(X_batch[i])
        if len(X_batch) < self.batch_size:
            diff = self.batch_size-len(X_batch)
            for _ in range(diff):
                X_batch.append(
                    np.zeros((self.target_size[0], self.target_size[1], 3), dtype=np.uint8))
        return np.array(X_batch, dtype=np.uint8)[:self.batch_size]


# %%
UnitTest = False
if UnitTest:
    img = HistoImage(
        "../Prostate/Slide-01.svs", verbose=1299)
    print(len(img))
    for batch in img:
        for i in range(len(batch)):
            plt.imshow(batch[i])
            plt.show()
    print(img.ThumbnailMask.shape)
    print(img.ThumbnailImage.size)
    print(img.store_coordination)
    print(img.dimension)

# %%
