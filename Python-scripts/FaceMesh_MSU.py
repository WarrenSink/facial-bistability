#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 11:28:47 2020

@author: Warren.Sink
"""
import sys
import os
import numpy as np
import pandas as pd
import cv2
import mediapipe as mp

mp_drawing = mp.solutions.drawing_utils
mp_face_mesh = mp.solutions.face_mesh

#get path to where this python file is contained
base_dir = os.path.dirname(__file__)

#get cli argument where relative path is from base_dir
image_folder = sys.argv[1]
annotated_image_folder = sys.argv[2]

if not os.path.exists(annotated_image_folder):
        
    print("New directory created")
        
    os.makedirs(annotated_image_folder)
    
#create an empty list
filenames = []

#create empty list for names to which coordinates belong
img_names = []

#create empty list for 
no_face_detected = []

#fill the empty file with the new filenames 
for idx, file in enumerate(os.listdir(base_dir + image_folder)):
    filename = os.fsdecode(file)
    if filename.endswith(".jpg"): 
        filenames.append(filename)
        continue
    else:
        continue    

#sort the filenames so they're in alphabetical order    
filenames.sort()

# script source: https://google.github.io/mediapipe/solutions/face_mesh
# For static images:
landmark_list = []
face_mesh = mp_face_mesh.FaceMesh(
    static_image_mode=True,
    max_num_faces=1,
    min_detection_confidence=0.5)
drawing_spec = mp_drawing.DrawingSpec(thickness=1, circle_radius=1)

#get results of dnn finding face and corresponding landmarks
#and turn results into a string in a text file
with open(annotated_image_folder + '.txt','w') as f:
    
  for idx, file in enumerate(filenames):
      
    print(base_dir + image_folder + '/' +file)
    image = cv2.imread(base_dir + image_folder + '/' + file)
    img_names.append(file)
  
    results = face_mesh.process(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))

    if not results.multi_face_landmarks:
      no_face_detected.append(file)
      continue
    annotated_image = image.copy()

    landmarks = results.multi_face_landmarks  
      
    for idy, face_landmarks in enumerate(results.multi_face_landmarks):
      f.write(str(face_landmarks))
      mp_drawing.draw_landmarks(
          image=annotated_image,
          landmark_list=face_landmarks,
          connections=mp_face_mesh.FACE_CONNECTIONS,
          landmark_drawing_spec=drawing_spec,
          connection_drawing_spec=drawing_spec)
    cv2.imwrite('./' + annotated_image_folder + '/' + file[:5] + '.png', annotated_image)
  f.close()
face_mesh.close()

#in order only keep images with results 
#img_names is used to assign each set of 468 landmarks as belonging to certain image
img_names = [item for item in img_names if item not in no_face_detected]


#when parsing the string of landmarks of all the images,
#the landmarks occur in blocks of 5 with the first line being 
#something like "landmark {", which I delimit by "landmark",
#the next three lines are the x, y, and z coordinates, and the
#fifth line is delimited by the closing right bracket
landmarkOrBrack = ['landmark', '}'] 
coordinate_list = ['x:','y:','z:']
img_names_index = 0 
with open(annotated_image_folder+".txt", "r+") as f, open(annotated_image_folder+'parsed.txt', 'w') as newf:
    idz = -1    
    for line in f:
        
        idz = idz + 1
        if idz ==  0 or idz == 2340:
    
            newf.write(img_names[img_names_index] + "\n") 
            img_names_index = img_names_index + 1
            
            idz = 0
            newf.write("landmark_" + str(int((idz/5)+1)) + "\n")
            
        else:
            
            if not any(landmarkOrBrack in line for landmarkOrBrack in landmarkOrBrack):
                for word in coordinate_list:
                    line = line.replace(word, "")
                newf.write(line.lstrip())
            else:
                if '}' in line:
                    continue
                else:
                    newf.write("landmark_" + str(int((idz/5)+1)) + "\n")
