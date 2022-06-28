# Tools - Convert lecture videos to slides
**Language:** Python <br>
**Softwares / packages:** cv2, skimage <br>
## Video information
![videoinfo](Lecture%20videos%20to%20slides/videoinfo.png)<br>
## Pre-processing working directory
```<Python>
import cv2
import os
import shutil
import winsound
import  numpy as np
from skimage.measure import compare_ssim
import av

rawFilesFolder = "Path" # Subject to change
videoFolder = rawFilesFolder + "\\"
framesFolder = rawFilesFolder + "\\Frames\\"
keyFramesFolder = rawFilesFolder + "\\KeyFrames\\"

rawFilesList = os.listdir(rawFilesFolder)
videosList = []
# Check if the folder contains non mp4 files
for i in range(len(rawFilesList)):
    if rawFilesList[i][-4:]== ".mp4":
        videosList.append(rawFilesList[i])
    else:
        print(rawFilesList[i])
print("There are %s videos. Do manual check. " %(len(videosList)) )
videosList


# Truncate long file names
videosListTrim = []
for i in range(len(videosList)):
    if len(videosList[i]) > 120:
        print("--->" + str(len(videosList[i])) + "<---")
        newFileName = videosList[i][:50] + " ... " + videosList[i][-50:]
        videosListTrim.append(newFileName)
        print(newFileName)
    else:
        print(len(videosList[i]))
        videosListTrim.append(videosList[i])
```
## Extracting a frame every 3 seconds from the video
```<Python>
def get_source_info_opencv(source_name): 
    try:
        cap = cv2.VideoCapture(source_name)
        width = cap.get(cv2.CAP_PROP_FRAME_WIDTH )
        height = cap.get(cv2.CAP_PROP_FRAME_HEIGHT)
        fps = cap.get(cv2.CAP_PROP_FPS)
        num_frames = cap.get(cv2.CAP_PROP_FRAME_COUNT)
        print("width:{} \nheight:{} \nfps:{} \nnum_frames:{}".format(width, height, fps, num_frames))
    except (OSError, TypeError, ValueError, KeyError, SyntaxError) as e:
        print("init_source:{} error. {}\n".format(source_name, str(e)))
    return fps, num_frames


for i in range(len(videosList)):
    print("Processing No.", str(i+1), "video.")
    currSourceName = videoFolder + videosList[i]
    currFramesPath = framesFolder + videosListTrim[i]
    os.makedirs(currFramesPath)
    
    fps, num_frames = get_source_info_opencv(currSourceName)
    fpsRound = round(fps)
    framesPerThreeSec = fpsRound * 3
    
    videoContainer = av.open(currSourceName)
    currProgess = 0
    preProgess = 0
    for frame in videoContainer.decode(video = 0):
        if frame.index % framesPerThreeSec == 0:
            frame.to_image().save(currFramesPath + "\\frame-%09d.jpg" % frame.index)
        currProgess = round(frame.index / num_frames, 2)
        if (currProgess * 100 % 10 == 0) & (currProgess != preProgess):
            print(currProgess * 100, "%")
            preProgess = currProgess
```
![videoprocess](Lecture%20videos%20to%20slides/videoprocess.png)<br>
![frames](Lecture%20videos%20to%20slides/frames.png)<br>
## Extracting key frames from frames
```<Python>
class CompareImage():

    def compare_image(self, path_image1, path_image2):

        imageA = cv2.imread(path_image1)
        imageB = cv2.imread(path_image2)

        grayA = cv2.cvtColor(imageA, cv2.COLOR_BGR2GRAY)
        grayB = cv2.cvtColor(imageB, cv2.COLOR_BGR2GRAY)

        (score, diff) = compare_ssim(grayA, grayB, full=True)
        print("SSIM: {}".format(score))
        return score

framesFolderList = os.listdir(framesFolder)

for i in range(len(framesFolderList)):
    print("Processing No.", str(i+1), "frames folder.")
    currFramesFolder = framesFolder + framesFolderList[i]
    currKeyFramesPath = keyFramesFolder + framesFolderList[i]
    os.makedirs(currKeyFramesPath)
    filesList = os.listdir(currFramesFolder)
    for j in range(len(filesList) - 1):
        print("Progress: " + str(j/len(filesList)*100) + "%")
        imageFile1 = currFramesFolder + "\\" + filesList[j]
        imageFile2 = currFramesFolder + "\\" + filesList[j+1]
        compare_image = CompareImage()
        curr_score = compare_image.compare_image(imageFile1, imageFile2)
        if curr_score < 0.94: 
            shutil.copy(imageFile1, currKeyFramesPath + "\\" + filesList[j])
    shutil.copy(currFramesFolder + "\\" + filesList[-1], currKeyFramesPath + "\\" + filesList[-1])
```
![imageprocess](Lecture%20videos%20to%20slides/imageprocess.png)<br>
![keyframes](Lecture%20videos%20to%20slides/keyframes.png)<br>