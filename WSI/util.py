from fileinput import FileInput
from matplotlib.image import thumbnail
import openslide
from tqdm.std import tqdm
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.patches import Rectangle
import cv2
from skimage.measure import label, regionprops
from skimage.color import label2rgb
from joblib import Parallel, delayed
from PIL import Image, ImageFilter
import random
import os
import tensorflow as tf
import pandas as pd
import pyvips
random.seed(1234)  
plt.figure(figsize=(20,20)) 
TestUnit = False

fold="/home/truebinding/coro/H&E (Karun Mutreja 4)/"
if TestUnit:
    filenames=os.listdir(fold)
def GetMaskWithReference(ref, size=(512,512) ,threshold=2, channel=0,strategy=0):
    if isinstance(ref, openslide.OpenSlide):
        ref_img=np.array(ref.get_thumbnail(size))
    else:
        ref_img= ref
    ref_img=ref_img[...,channel]
    mask=np.zeros_like(ref_img)
    
    mask[ref_img>=threshold]=1
    if strategy==1:
        mask=cv2.dilate(mask, np.ones((2,2), np.uint8))
        mask[mask>=0.5]=1
        mask[mask<0.5]=0
    if strategy==2:
        mask=cv2.dilate(mask, np.ones((4,4), np.uint8))
        mask[mask>=0.5]=1
        mask[mask<0.5]=0
    if strategy==0:
        mask=np.array(Image.fromarray(mask).filter(ImageFilter.GaussianBlur(4)))
        mask[mask>=0.5]=1
        mask[mask<0.5]=0
    return mask
def GetMask(img,DebugMode=False):
    #img = openslide.OpenSlide(filename=filename)
    img = img.associated_images["thumbnail"]
    img_x=img.filter(ImageFilter.GaussianBlur(8))
    rgb = np.array(img_x.convert("RGB")).astype(np.uint8)
    hls = cv2.cvtColor(rgb, cv2.COLOR_RGB2HLS)
    mask = 255-cv2.inRange(hls, np.array([0,225,0]), np.array([255,255,255]))
    rgb = cv2.bitwise_and(rgb,rgb, mask= mask)
    #BLUE
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    lwr = np.array([90,50,50])
    upr = np.array([130,255,255])
    msk = cv2.inRange(hsv, lwr, upr)
    
    mask=mask-msk
    mask[mask<0]=0
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    hsv = cv2.bitwise_and(hsv,hsv, mask= mask)
    msk = cv2.inRange(hsv, (36, 0, 0), (90, 255,255))
    
    mask=mask-msk
    mask[mask<0]=0


    kernel=np.ones((5,5), np.uint8)
    mask=cv2.erode(mask, kernel, iterations=1)
    mask=cv2.dilate(mask,kernel, iterations=1)
    #img=np.array(img.convert("RGB"))
    #res = cv2.bitwise_and(img,img, mask= mask)
    if DebugMode:
        plt.imshow(mask)
        plt.show()
    return mask
def GetImagePatchCoordinationFromMaskFile(filename, 
                            MaskFilename,
                            MinFilledSize=250,
                            target_magnification=20,
                            ImageMagnification=10,
                            threshold=2,
                            target_patch_size=(512,512),
                            overlap_rate = (0.5,  0.5),
                            background_threshold=0.20,
                            background_color=200,
                            DebugMode=False,
                            DebugModeV2=False,
                            Limit=200, strategy=0):
    img=openslide.OpenSlide(filename)
    thumbnail_img=img.get_thumbnail((512,512))
    mask=GetMaskWithReference(openslide.OpenSlide(MaskFilename), threshold=threshold, strategy=strategy)
    zoom=ImageMagnification
    w,h=img.level_dimensions[0]
    img.close()
    factor = target_magnification/zoom
    thumbnail_img=np.array(thumbnail_img)
    h_tiny,w_tiny, _= thumbnail_img.shape
    #print(h_tiny,w_tiny)
    label_image = label(mask)
    regions_bb= []
    for region in regionprops(label_image):
        # take regions with large enough areas
        if region.filled_area >= MinFilledSize:
            regions_bb.append(region)
    if DebugMode:
        plt.imshow(thumbnail_img)
        plt.show()
    store_coordination = []
    for region in regions_bb:
        reg=region.bbox
        min_row, min_col, max_row,  max_col=reg 
        row_l_ = max_row - min_row
        col_l_ = max_col - min_col
        
        #x,y coordinate top_left
        Pt_top_left=((min_col/w_tiny)*w*factor,(min_row/h_tiny)*h*factor)
        #x,y coordinate bottom right
        Pt_bottom_right = ((max_col/w_tiny)*w*factor,(max_row/h_tiny)*h*factor)
        org_x_l_ = Pt_bottom_right[0]-Pt_top_left[0]
        org_y_l_ =  Pt_bottom_right[1]-Pt_top_left[1]
        number_of_patches_x_axis=org_x_l_/target_patch_size[0]
        number_of_patches_y_axis=org_y_l_/target_patch_size[1]
        
        patch_size_in_tiny_y_axis=round(row_l_/number_of_patches_y_axis)
        patch_size_in_tiny_x_axis=round(col_l_/number_of_patches_x_axis)
        
        if patch_size_in_tiny_y_axis!=patch_size_in_tiny_x_axis:
            if patch_size_in_tiny_y_axis<patch_size_in_tiny_x_axis:
                patch_size_in_tiny_x_axis=patch_size_in_tiny_y_axis
            else:
                patch_size_in_tiny_y_axis=patch_size_in_tiny_x_axis

        masks = np.zeros((thumbnail_img.shape[0],thumbnail_img.shape[1]))

        masks[min_row:max_row,min_col:max_col]=region.image

        thumbnail_img_corrected = thumbnail_img.copy()
        for i in range(3):
            thumbnail_img_corrected[...,i]=thumbnail_img_corrected[...,i]*(masks)
        for y in range(min_row,max_row-math.floor(patch_size_in_tiny_y_axis*overlap_rate[1]), math.floor(patch_size_in_tiny_y_axis*overlap_rate[1])):
            for x in range(min_col, max_col-math.floor(patch_size_in_tiny_y_axis*overlap_rate[0]), math.floor(patch_size_in_tiny_x_axis*overlap_rate[0])):
                
                ic= thumbnail_img_corrected[y:y+patch_size_in_tiny_y_axis, x:x+patch_size_in_tiny_x_axis]
                img_hsv=cv2.cvtColor(ic, cv2.COLOR_RGB2HSV)
                bg_img=cv2.cvtColor(ic,cv2.COLOR_RGB2GRAY)

                total_count = patch_size_in_tiny_x_axis*patch_size_in_tiny_y_axis
                # lower mask (0-10)
                lower_red = np.array([0,50,50])
                upper_red = np.array([10,255,255])
                mask0 = cv2.inRange(img_hsv, lower_red, upper_red)

                # upper mask (170-180)
                lower_red = np.array([170,50,50])
                upper_red = np.array([180,255,255])
                mask1 = cv2.inRange(img_hsv, lower_red, upper_red)
                mask_b = mask0+mask1
                
                mask_b[mask_b>0]=1
                non_z_count_b = np.count_nonzero(mask_b)
                per_red_background = non_z_count_b/total_count
                if per_red_background>0.5:

                    continue

                if DebugModeV2:
                    plt.imshow(bg_img, cmap="gray")
                    plt.show()
                    plt.imshow(bg_img>background_color, cmap="gray")
                    plt.show()
                non_z_count_ = np.count_nonzero(bg_img==0)
                non_z_count = np.count_nonzero(bg_img>background_color)+non_z_count_
                per_white_background = non_z_count/total_count
                
                #if DebugMode and per_white_background<background_threshold:
                #    ax.add_patch(Rectangle((x,y),patch_size_in_tiny_y_axis, patch_size_in_tiny_x_axis, edgecolor='r', facecolor="none"))
                if per_white_background<background_threshold:
                    y_org_location = int(round(y/h_tiny*h))
                    x_org_location = int(round(x/w_tiny*w))
                    store_coordination.append([x_org_location,y_org_location,target_patch_size[0]*(zoom//target_magnification),target_patch_size[1]*(zoom//target_magnification) ])
    if DebugMode:
            fig,ax = plt.subplots()
            ax.imshow(thumbnail_img)
            
            for coord in store_coordination:
                x=round((coord[0]/w)*w_tiny)
                y=round((coord[1]/h)*h_tiny)
                ax.add_patch(Rectangle((x,y),patch_size_in_tiny_y_axis, patch_size_in_tiny_x_axis, edgecolor='r', facecolor="none"))
            plt.show()
    if Limit!=0:
        if len(store_coordination)>Limit:
            store_coordination=random.sample(store_coordination,Limit)
    return store_coordination


def GetImagePatchCoordination(filename, 
                            MinFilledSize=250,
                            target_magnification=20,
                            target_patch_size=(512,512),
                            overlap_rate = (0.99,  0.99),
                            background_threshold=0.75,
                            background_color=220,
                            DebugMode=False,
                            DebugModeV2=False,
                            Limit=200
):
    img=openslide.OpenSlide(filename)
    thumbnail_img=img.associated_images["thumbnail"].copy()
    mask=GetMask(img,DebugMode)
    zoom=int(img.properties["aperio.AppMag"])
    #print(zoom)
    w,h=img.level_dimensions[0]
    #print(w,h)
    img.close()
    factor = target_magnification/zoom
    #image.crop(x, y, size, size)
    #print(factor)
    thumbnail_img=np.array(thumbnail_img)
    h_tiny,w_tiny, _= thumbnail_img.shape
    #print(h_tiny,w_tiny)
    label_image = label(mask)
    
    regions_bb= []
    for region in regionprops(label_image):
        # take regions with large enough areas
        if region.filled_area >= MinFilledSize:
            regions_bb.append(region)
    if DebugMode:
        plt.imshow(thumbnail_img)
        plt.show()
    #print(len(regions_bb))
    
    store_coordination = []
    for region in regions_bb:
        reg=region.bbox
        min_row, min_col, max_row,  max_col=reg 
        row_l_ = max_row - min_row
        col_l_ = max_col - min_col
        
        #x,y coordinate top_left
        Pt_top_left=((min_col/w_tiny)*w*factor,(min_row/h_tiny)*h*factor)
        #x,y coordinate bottom right
        Pt_bottom_right = ((max_col/w_tiny)*w*factor,(max_row/h_tiny)*h*factor)
        org_x_l_ = Pt_bottom_right[0]-Pt_top_left[0]
        org_y_l_ =  Pt_bottom_right[1]-Pt_top_left[1]
        number_of_patches_x_axis=org_x_l_/target_patch_size[0]
        number_of_patches_y_axis=org_y_l_/target_patch_size[1]
        
        patch_size_in_tiny_y_axis=round(row_l_/number_of_patches_y_axis)
        patch_size_in_tiny_x_axis=round(col_l_/number_of_patches_x_axis)
        
        if patch_size_in_tiny_y_axis!=patch_size_in_tiny_x_axis:
            if patch_size_in_tiny_y_axis<patch_size_in_tiny_x_axis:
                patch_size_in_tiny_x_axis=patch_size_in_tiny_y_axis
            else:
                patch_size_in_tiny_y_axis=patch_size_in_tiny_x_axis

        masks = np.zeros((thumbnail_img.shape[0],thumbnail_img.shape[1]))
        #print(masks.shape)
        #print(region.image.shape,row_l_,col_l_)
        
        masks[min_row:max_row,min_col:max_col]=region.image
        #plt.imshow(masks)
        #plt.show()
        thumbnail_img_corrected = thumbnail_img.copy()
        for i in range(3):
            thumbnail_img_corrected[...,i]=thumbnail_img_corrected[...,i]*(masks)
        for y in range(min_row,max_row-math.floor(patch_size_in_tiny_y_axis*overlap_rate[1]), math.floor(patch_size_in_tiny_y_axis*overlap_rate[1])):
            for x in range(min_col, max_col-math.floor(patch_size_in_tiny_y_axis*overlap_rate[0]), math.floor(patch_size_in_tiny_x_axis*overlap_rate[0])):
                
                ic= thumbnail_img_corrected[y:y+patch_size_in_tiny_y_axis, x:x+patch_size_in_tiny_x_axis]
                img_hsv=cv2.cvtColor(ic, cv2.COLOR_RGB2HSV)
                bg_img=cv2.cvtColor(ic,cv2.COLOR_RGB2GRAY)

                total_count = patch_size_in_tiny_x_axis*patch_size_in_tiny_y_axis
                # lower mask (0-10)
                lower_red = np.array([0,50,50])
                upper_red = np.array([10,255,255])
                mask0 = cv2.inRange(img_hsv, lower_red, upper_red)

                # upper mask (170-180)
                lower_red = np.array([170,50,50])
                upper_red = np.array([180,255,255])
                mask1 = cv2.inRange(img_hsv, lower_red, upper_red)
                mask_b = mask0+mask1
                
                mask_b[mask_b>0]=1
                non_z_count_b = np.count_nonzero(mask_b)
                per_red_background = non_z_count_b/total_count
                #plt.imshow(mask_b)
                #plt.show()
                if per_red_background>0.5:
                    #plt.imshow(img_hsv[...,0])
                    #plt.show()
                    #print(non_z_count_b)
                    continue

                if DebugModeV2:
                    plt.imshow(bg_img, cmap="gray")
                    plt.show()
                    plt.imshow(bg_img>background_color, cmap="gray")
                    plt.show()
                non_z_count_ = np.count_nonzero(bg_img==0)
                non_z_count = np.count_nonzero(bg_img>background_color)+non_z_count_
                per_white_background = non_z_count/total_count
                
                #if DebugMode and per_white_background<background_threshold:
                #    ax.add_patch(Rectangle((x,y),patch_size_in_tiny_y_axis, patch_size_in_tiny_x_axis, edgecolor='r', facecolor="none"))
                if per_white_background<background_threshold:
                    y_org_location = int(round(y/h_tiny*h))
                    x_org_location = int(round(x/w_tiny*w))
                    store_coordination.append([x_org_location,y_org_location,target_patch_size[0]*(zoom//target_magnification),target_patch_size[1]*(zoom//target_magnification) ])
    if DebugMode:
            fig,ax = plt.subplots()
            ax.imshow(thumbnail_img)
            
            for coord in store_coordination:
                x=round((coord[0]/w)*w_tiny)
                y=round((coord[1]/h)*h_tiny)
                ax.add_patch(Rectangle((x,y),patch_size_in_tiny_y_axis, patch_size_in_tiny_x_axis, edgecolor='r', facecolor="none"))
            plt.show()
    if Limit!=0:
        if len(store_coordination)>Limit:
            store_coordination=random.sample(store_coordination,Limit)
    return store_coordination
#%%
if TestUnit:
    from tqdm import tqdm
    for fl_ in tqdm(filenames[50:60]):
        filename=f"{fold}{fl_}"
        print(len(GetImagePatchCoordination(filename=filename, DebugMode=TestUnit)))
# %%
class OpenSlide(tf.keras.utils.Sequence):
    def __init__(self, csv_filename : str, column_source : str, endpoint : str,batch_size : int, target_size : list, target_magnification:int, augmentation):
        '''
        @par csv_filename: the csv file should include columns for filename, endpoint, x,y.
        @par column_source
        @par endpoint
        @par batch_size
        @par target_size
        @par augmentation
        '''
        self.csv_filename=csv_filename
        self.column_source=column_source
        self.endpoint=endpoint
        self.batch_size=batch_size
        self.regions = []
        self.data = pd.read_csv(csv_filename)
        self.filenames=self.data[self.column_source].tolist()
        self.outcome = self.data[self.endpoint].tolist()
        self.X_point = self.data['x'].tolist()
        self.Y_point = self.data['y'].tolist()
        self.target_size=target_size
        self.augmentation=augmentation
        self.target_magnification=target_magnification
        self.open_slides = {}
    def __len__(self):
         return math.ceil(len(self.filenames) / self.batch_size)
    def GetImage(self, filename, region):
        format_to_dtype = {
            'uchar': np.uint8,
            'char': np.int8,
            'ushort': np.uint16,
            'short': np.int16,
            'uint': np.uint32,
            'int': np.int32,
            'float': np.float32,
            'double': np.float64,
            'complex': np.complex64,
            'dpcomplex': np.complex128,
        }
        image = pyvips.Image.new_from_file(filename)
        img=openslide.OpenSlide(filename)
        zoom=int(img.properties["aperio.AppMag"])
        img.close()
        size=self.target_size[0]*(zoom//self.target_magnification)
        img=image.crop(region[0], region[1], size, size)#location=region, level=0, size=(size,size))
        img = np.ndarray(buffer=img.write_to_memory(),
                   dtype=format_to_dtype[img.format],
                   shape=[img.height, img.width, img.bands])
        img=Image.fromarray(img)
        img=img.convert("RGB").resize(self.target_size)
        return np.array(img)

    def __getitem__(self, idx):
        batch_filenames = self.filenames[idx * self.batch_size:(idx + 1) *
        self.batch_size]

        batch_X_point = self.X_point[idx * self.batch_size:(idx + 1) *
        self.batch_size]

        batch_Y_point = self.Y_point[idx * self.batch_size:(idx + 1) *
        self.batch_size]

        batch_y = self.outcome[idx * self.batch_size:(idx + 1) *
        self.batch_size]
        
        batch_x=[]
        for i in range(len(batch_filenames)):
            filename=batch_filenames[i]
            region=(batch_X_point[i], batch_Y_point[i])
            img=self.GetImage(filename,region)
            batch_x.append(img)

        batch_x=np.array(batch_x)

        if self.augmentation is not None:
            batch_x = self.augmentation(images=batch_x)

        return batch_x, np.array(batch_y)
        
class ImageGeneratorCSV(tf.keras.utils.Sequence):
    def __init__(self, csv_filename : str, source : str, filename : str, endpoint : str,batch_size : int, target_size : list, target_magnification:int, augmentation):
        '''
        @par csv_filename: the csv file should include columns for filename, endpoint, x,y.
        @par column_source
        @par endpoint
        @par batch_size
        @par target_size
        @par augmentation
        '''
        self.csv_filename=csv_filename
        self.endpoint=endpoint
        self.batch_size=batch_size
        self.regions = []
        self.data = pd.read_csv(csv_filename)
        self.filenames=self.data[filename].tolist()
        self.outcome = self.data[endpoint].tolist()
        self.target_size=target_size
        self.augmentation=augmentation
        self.source=source
        self.n_class = len(np.unique(self.outcome))
    def __len__(self):
         return math.floor(len(self.filenames) / self.batch_size)
    def GetImage(self, filename):
        img=Image.open(f"{self.source}/{filename}")
        if img.size!=self.target_size:
            img=img.resize(self.target_size)

        return np.array(img)

    def __getitem__(self, idx):
        batch_filenames = self.filenames[idx * self.batch_size:(idx + 1) *
        self.batch_size]
        batch_y = self.outcome[idx * self.batch_size:(idx + 1) *
        self.batch_size]
        batch_x=[]
        for i in range(len(batch_filenames)):
            filename=batch_filenames[i]
            img=self.GetImage(filename)
            #print(img.shape,filename)
            '''
            if img is None:
                img = np.zeros((self.target_size,self.target_size,3))
            if img.shape[0]==0:
                img = np.zeros((self.target_size,self.target_size,3))
            if img.shape[1]==0:
                img = np.zeros((self.target_size,self.target_size,3))
            '''
            batch_x.append(img)

        batch_x=np.array(batch_x)

        if self.augmentation is not None:
            batch_x = self.augmentation(images=batch_x)

        return batch_x, tf.keras.utils.to_categorical(np.array(batch_y), self.n_class)