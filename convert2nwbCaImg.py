'''
Convert two-photon calcium imaging data to the NWB format in Python

Run this script to convert two-photon calcium imaging data and associated
intracellular physiology data generated at the University of Bristol
(UoB) to the Neurodata Without Borders (NWB) file format. This script is
explained in the accompanying Bristol GIN for Calcium Imaging Data
tutorial available at
https://dervinism.github.io/bristol-neuroscience-data-guide/tutorials/Bristol%20GIN%20for%20Calcium%20Imaging%20Data.html

You can use this script to get an idea of how to convert your own
optical imaging data to the NWB file format.
'''

import scipy.io
import numpy as np
from datetime import datetime

from pynwb import NWBFile, NWBHDF5IO
from pynwb.core import DynamicTable
from pynwb.file import Subject
from pynwb.ophys import OpticalChannel
from pynwb.base import Images
from pynwb.image import RGBImage, GrayscaleImage

from localFunctions import appendLinescans, reshapeData, setTwoPhotonSeries, setDeltaFSeries, setCClampSeries

# Record metadata
# Project (experiment) metadata:
projectName = 'Intracellular Ca2+ dynamics during plateau potentials trigerred by Schaffer collateral stimulation'
experimenter = 'Matt Udakis'
institution = 'University of Bristol'
publications = 'In preparation'
lab = 'Jack Mellor lab'
brainArea = 'Hippocampus CA1-2'
greenIndicator = 'Fluo5f'
redIndicator = 'Alexa594'

# Animal metadata
animalID = 'm1'
ageInDays = 100
age = 'P'+str(ageInDays)+'D' # Convert to ISO8601 format: https://en.wikipedia.org/wiki/ISO_8601#Durations
strain = 'C57BL/6J'
sex = 'M'
species = 'Mus musculus'
weight = []
description = '001' # Animal testing order.

# Session metadata
startYear = 2020
startMonth = 12
startDay = 4
startTime = datetime(startYear, startMonth, startDay)
year = str(startYear); year = year[2:4]
month = str(startMonth)
if len(month) == 1:
  month = '0'+month
day = str(startDay)
if len(day) == 1:
  day = '0'+day
sliceNumber = 2
cellNumber = 1
imagingRate = 1/21 # A single linescan duration is 1sec with 20sec period in between linescans
lineRate = 1000. # A number of lines scanned in a second.
sessionID = animalID + '_' + year + month + day + '_s' + str(sliceNumber) + '_c' + str(cellNumber) # mouse-id_time_slice-id_cell-id
sessionDescription = 'Single cell imaging in a slice combined with somatic current clamp recordings and the stimulation of Schaffer collaterals'
sessionNotes = '201204 - Slice 2 Imaging 3 dendrite regions roughly in the SO SR and SLM regions   ' +\
               'Same line scan with same intensity of stim (2.3V) at different locations along the cell   ' +\
               'Ephys frames to match up with linescans   ' +\
               'Frames   ' +\
               '1-8 Bottom Den   ' +\
               '10-19 Middle Den   ' +\
               '22-28 Top Dendrite   ' +\
               'Missed a few of the imaging with the ephys so more Ephys traces than linescans   ' +\
               'by the end of the experiment the top or neuron started to bleb.'

# Assign NWB file fields
nwb = NWBFile(
  session_description = sessionDescription,
  identifier = sessionID, 
  session_start_time = startTime, 
  experimenter = experimenter,  # optional
  session_id = sessionID,  # optional
  institution = institution,  # optional
  related_publications = publications,  # optional
  notes = sessionNotes,  # optional
  lab = lab) # optional

# Create subject object
nwb.subject = Subject(
  subject_id = animalID,
  age = age,
  description = description,
  species = species,
  sex = sex)

# Convert calcium imaging data
# Create optical channels
green_optical_channel = OpticalChannel(
  name = "OpticalChannel",
  description = 'green channel corresponding to ' + greenIndicator,
  emission_lambda=516.)

red_optical_channel = OpticalChannel(
  name = "OpticalChannel",
  description = 'red channel corresponding to ' + redIndicator,
  emission_lambda=616.)

# Create the imaging device object
device = nwb.create_device(
    name = '2P_microscope',
    description = 'Two-photon microscope',
    manufacturer = 'Scientifica')

# Create imaging plane objects
green_imaging_plane = nwb.create_imaging_plane(
  name = 'green_imaging_plane',
  optical_channel = green_optical_channel,
  imaging_rate = imagingRate,
  description = 'The plane for imaging calcium indicator Fluo5f.',
  device = device,
  excitation_lambda = 810.,
  indicator = 'Fluo5f',
  location = brainArea)

red_imaging_plane = nwb.create_imaging_plane(
  name = 'red_imaging_plane',
  optical_channel = red_optical_channel,
  imaging_rate = imagingRate,
  description = 'The plane for imaging calcium indicator Alexa594.',
  device = device,
  excitation_lambda = 810.,
  indicator = 'Alexa594',
  location = brainArea)

# Load data
dendrite = 'Bot'
data1 = scipy.io.loadmat('../Analysed/' + year + month + day + '__s' + str(sliceNumber) + \
  'd' + str(cellNumber) + '_004_ED__1 ' + dendrite + 'den_Analysed.mat', squeeze_me=True)['Analysed_data']
dendrite = 'Mid'
data2 = scipy.io.loadmat('../Analysed/' + year + month + day + '__s' + str(sliceNumber) + \
  'd' + str(cellNumber) + '_004_ED__1 ' + dendrite + 'den_Analysed.mat', squeeze_me=True)['Analysed_data']
dendrite = 'Top'
data3 = scipy.io.loadmat('../Analysed/' + year + month + day + '__s' + str(sliceNumber) + \
  'd' + str(cellNumber) + '_004_ED__1 ' + dendrite + 'den_Analysed.mat', squeeze_me=True)['Analysed_data']

# Append raw linescans so they have equal widths
data1['Flur5_denoised'].fill(appendLinescans(data1['Flur5_denoised'].item()))
data1['Alexa_denoised'].fill(appendLinescans(data1['Alexa_denoised'].item()))
data2['Flur5_denoised'].fill(appendLinescans(data2['Flur5_denoised'].item()))
data2['Alexa_denoised'].fill(appendLinescans(data2['Alexa_denoised'].item()))
data3['Flur5_denoised'].fill(appendLinescans(data3['Flur5_denoised'].item()))
data3['Alexa_denoised'].fill(appendLinescans(data3['Alexa_denoised'].item()))
nFrames = [len(data1['Flur5_denoised'].item()), \
  len(data2['Flur5_denoised'].item()), len(data1['Flur5_denoised'].item())]

# Add optical physiology data: Fluorescence
input = {
  'indicator': greenIndicator,
  'imagingPlane': green_imaging_plane,
  'imagingRate': imagingRate,
  'lineRate': lineRate,
  'data': reshapeData(data1['Flur5_denoised'].item()),
  'dendriteID': 'bottom'}
nwb = setTwoPhotonSeries(nwb, input)

input['data'] = reshapeData(data2['Flur5_denoised'].item())
input['dendriteID'] = 'middle'
nwb = setTwoPhotonSeries(nwb, input)

input['data'] = reshapeData(data3['Flur5_denoised'].item())
input['dendriteID'] = 'top'
nwb = setTwoPhotonSeries(nwb, input)

input['indicator'] = redIndicator
input['imagingPlane'] = red_imaging_plane
input['data'] = reshapeData(data1['Alexa_denoised'].item())
input['dendriteID'] = 'bottom'
nwb = setTwoPhotonSeries(nwb, input)

input['data'] = reshapeData(data2['Alexa_denoised'].item())
input['dendriteID'] = 'middle'
nwb = setTwoPhotonSeries(nwb, input)

input['data'] = reshapeData(data3['Alexa_denoised'].item())
input['dendriteID'] = 'top'
nwb = setTwoPhotonSeries(nwb, input)

# Add optical physiology data: Delta fluorescence
input = {
  'indicator': redIndicator,
  'imagingPlane': red_imaging_plane,
  'imagingRate': imagingRate,
  'lineRate': lineRate,
  'data': np.expand_dims(np.transpose(data1['Calcium_deltaF'].item()), axis=2),
  'dendriteID': 'bottom'}
nwb = setDeltaFSeries(nwb, input)

input['data'] = np.expand_dims(np.transpose(data2['Calcium_deltaF'].item()), axis=2),
input['dendriteID'] = 'middle'
nwb = setDeltaFSeries(nwb, input)

input['data'] = np.expand_dims(np.transpose(data3['Calcium_deltaF'].item()), axis=2),
input['dendriteID'] = 'top'
nwb = setDeltaFSeries(nwb, input)

# Add images
neuron_image = RGBImage(
  name = 'neuron_image',
  data = data1['Neuron'].item(), # required: [height, width, colour]
  description = 'RGB image of the full neuron.')

bottom_dend_image = GrayscaleImage(
  name = 'dendrite1_image',
  data = data1['ROI_img'].item(), # required: [height, width]
  description = 'Grayscale image of the bottom dendrite.')

middle_dend_image = GrayscaleImage(
  name = 'dendrite2_image',
  data = data2['ROI_img'].item(), # required: [height, width]
  description = 'Grayscale image of the middle dendrite.')

top_dend_image = GrayscaleImage(
  name = 'dendrite3_image',
  data = data3['ROI_img'].item(), # required: [height, width]
  description = 'Grayscale image of the top dendrite.')

image_collection = Images(
  name = 'ImageCollection',
  images = [neuron_image, bottom_dend_image, middle_dend_image, top_dend_image],
  description = 'A collection of neuron and dendrite images.')

nwb.add_acquisition(image_collection)

# Convert intracellular electrophysiology data
# Create the recording device object
device = nwb.create_device(
  name = 'Amplifier_Multiclamp_700A',
  description = 'Amplifier for recording current clamp data.',
  manufacturer = 'Molecular Devices')

electrode = nwb.create_icephys_electrode(
  name = 'icephys_electrode',
  description = 'A patch clamp electrode',
  location = 'Cell soma in CA1-2 of hippocampus',
  slice = 'slice #' + str(sliceNumber),
  device = device)

# Add current clamp data
input = {
  'ephysTime': data1['Ephys_Time'],
  'nSweeps': nFrames,
  'data': np.transpose(data1['Ephys_data'].item()),
  'imagingRate': imagingRate,
  'electrode': electrode,
  'dendriteID': 'bottom'}
nwb = setCClampSeries(nwb, input)

input['ephysTime'] = data2['Ephys_Time']
input['data'] = np.transpose(data2['Ephys_data'].item())
input['dendriteID'] = 'middle'
nwb = setCClampSeries(nwb, input)

input['ephysTime'] = data3['Ephys_Time']
input['data'] = np.transpose(data3['Ephys_data'].item())
input['dendriteID'] = 'top'
nwb = setCClampSeries(nwb, input)

# Save the converted NWB file
with NWBHDF5IO(sessionID + '.nwb', "w") as io:
  io.write(nwb)