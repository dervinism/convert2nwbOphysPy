'''
Local functions for convert2nwbCaImg
'''

import numpy as np
from pynwb.core import VectorData
from pynwb.ophys import TwoPhotonSeries
from pynwb.icephys import CurrentClampSeries


def appendLinescans(linescans):
  '''
  appendedLinescans = appendLinescans(linescans)
  
  Function takes linescans of unequal widths and appends them to have the
  the same width and outputs the appended linescans as a 3D matrix with
  the first dimension corresponding to the linescan number.
  '''

  # Find the largest width
  nScans = len(linescans)
  length = len(linescans[0])
  maxWidth = len(linescans[0][0])
  for iScan in range(nScans):
    maxWidth = max([maxWidth, len(linescans[iScan][0])])

  for iScan in range(nScans):
    array2add = np.empty([len(linescans[iScan]),maxWidth - len(linescans[iScan][0])])
    array2add[:] = np.nan
    linescans[iScan] = np.append(linescans[iScan], array2add, axis=1)
  
  return linescans


def setTwoPhotonSeries(nwb, input):
  '''
  nwb = setTwoPhotonSeries(nwb, input)
  
  Function creates, names, and adds a two-photon series data to a given NWB
  file.
  Input: nwb - the NWB file object.
         input - a dictionary with the following fields:
           indicator - a string variable with the calcium indicator name
                       (e.g., 'Fluo5f' or 'Alexa594').
           imagingPlane - the imaging plane object corresponding to the
                          indicator.
           imagingRate - a scalar with the rate at which individual full
                         linescans are obtained.
           lineRate - a scalar wit the frequency of individual lines within
                      a single linescan.
           data - an n-dimensional matrix containing two-photon series
                  linescan data. The first dimension corresponds to time
                  (that is, individual linescans). The second dimension
                  corresponds to individual lines along the vertical
                  denditic segment. The third dimension corresponds to the
                  dendritic width. NaNs are appended at the end of linescans
                  to make them equal in width.
           dendriteID - a string variable with the dendritic ID (i.e.,
                        'bottom', 'middle', and 'top').
  Output: nwb - the NWB file object containing the newly added two-photon
                series data.
  '''
  
  # Name the two-photon series to the NWB file
  if input['indicator'] in 'Fluo5f':
    opticalChannel = 'Green'
  elif input['indicator'] in 'Alexa594':
    opticalChannel = 'Red'
  
  if input['dendriteID'] in 'bottom':
    dendriteID = '1'
  elif input['dendriteID'] in 'middle':
    dendriteID = '2'
  elif input['dendriteID'] in 'top':
    dendriteID = '3'

  # Create image series
  image_series = TwoPhotonSeries(
    name = 'TwoPhotonSeries' + opticalChannel + dendriteID,
    description = input['indicator'] + ' linescans of the ' + input['dendriteID'] + ' dendrite',
    imaging_plane = input['imagingPlane'],
    starting_time = 0.0,
    rate = input['imagingRate'],
    scan_line_rate = input['lineRate'],
    data = input['data'],
    unit =  'a.u.',
    comments = 'This two-photon series contains ' + input['indicator'] + ' linescans of the ' +\
               input['dendriteID'] + ' (ROI) with the first dimension corresponding to time ' +\
               '(or to individual linescans). Each linescan is 1-sec in duration with ' +\
               '20-sec intervals between two linescans. The second dimension corresponds ' +\
               'to individual lines spanning the length of the dendrite in the ROI. ' +\
               'The third dimension corresponds to the width of the dendrite. ' +\
               'Some linescans may contain appended NaN values to make ' +\
               'widths of different linescans be equal within the same ROI.   ' +\
               'data_continuity = step')

  nwb.add_acquisition(image_series)
  return nwb


def setDeltaFSeries(nwb, input):
  '''
  nwb = setDeltaFSeries(nwb, input)
  
  Function creates, names, and adds a two-photon series delta fluorescence
  data to a given NWB file.
  Input: nwb - the NWB file object.
         input - a dictionary with the following fields:
           indicator - a string variable with the calcium indicator name
                       (e.g., 'Fluo5f' or 'Alexa594').
           imagingPlane - the imaging plane object corresponding to the
                          indicator.
           imagingRate - a scalar with the rate at which individual full
                         linescans are obtained.
           lineRate - a scalar wit the frequency of individual lines within
                      a single linescan.
           data - a 2D matrix containing two-photon series delta F data.
                  The first dimension corresponds to time (that is, individual
                  linescans). The second dimension corresponds to individual
                  lines along the vertical dendritic segment. The data is
                  averaged across the dendritic width.
           dendriteID - a string variable with the dendritic ID (i.e.,
                        'bottom', 'middle', and 'top').
  Output: nwb - the NWB file object containing the newly added two-photon
                delta F series data.
  '''

  # Name and add the two-photon delta F series to the NWB file
  if input['dendriteID'] in 'bottom':
    dendriteID = '1'
  elif input['dendriteID'] in 'middle':
    dendriteID = '2'
  elif input['dendriteID'] in 'top':
    dendriteID = '3'
  
  # Create image series
  image_series = TwoPhotonSeries(
    name = 'TwoPhotonDeltaFSeries' + dendriteID,
    description = 'Delta F data for the ' + input['dendriteID'] +\
                  ' calculated based on ' + input['indicator'] + '.',
    imaging_plane = input['imagingPlane'],
    starting_time = 0.0,
    rate = input['imagingRate'],
    scan_line_rate = input['lineRate'],
    data = input['data'],
    unit = 'normalised',
    comments = 'This two-photon series contains delta F data calculated based on ' +\
               input['indicator'] + ' for the ' + input['dendriteID'] + ' (ROI) ' +\
               'with the first dimension corresponding to time ' +\
               '(or to individual linescans). Each linescan is 1-sec in duration with ' +\
               '20-sec intervals between two linescans. The second dimension corresponds ' +\
               'to individual lines spanning the length of the dendrite in the ROI. ' +\
               'The data is averaged across the dendritic width.   ' +\
               'data_continuity = step')

  nwb.add_acquisition(image_series)
  return nwb


def setCClampSeries(nwb, input):
  '''
  nwb = setCClampSeries(nwb, input)
  
  Function creates, names, and adds a current clamp series data to a given
  NWB file.
  Input: nwb - the NWB file object.
         input - a dictionary with the following fields:
           ephysTime - a time vector for a single recording sweep.
           nSweeps - The total number of sweeps for every ROI.
           imagingRate - a scalar with the rate at which individual full
                         linescans are obtained.
           data - a 2D matrix containing somatic current clamp recordings
                  corresponding to the initial part of the calcium
                  imaging period at dendritic ROI. The first dimension
                  corresponds to individual recording sweeps. The second
                  dimension corresponds to individual sweep data samples.
           electrode - an electrode object.
           dendriteID - a string variable with the dendritic ID (i.e.,
                        'bottom', 'middle', and 'top').
  Output: nwb - the NWB file object containing the newly added current
                clamp series data.
  '''
  
  if input['dendriteID'] in 'bottom':
    nSweeps = input['nSweeps'][0]
    sweeps = list(range(0, nSweeps))
    dendriteID = '1'
  elif input['dendriteID'] in 'middle':
    nSweeps = input['nSweeps'][1]
    sweeps = list(range(input['nSweeps'][0],input['nSweeps'][0]+nSweeps))
    dendriteID = '2'
  else:
    nSweeps = input['nSweeps'][2]
    sweeps = list(range(input['nSweeps'][0]+input['nSweeps'][1],input['nSweeps'][0]+input['nSweeps'][1]+nSweeps))
    dendriteID = '3'

  timestamps = input['ephysTime']/1000; # sec

  current_clamp_series = CurrentClampSeries(
    name = 'CurrentClampSeries' + dendriteID,
    description = 'Somatic current clamp recording corresponding to ' +\
                  'the initial part of the calcium imaging period ' +\
                  'at ' + input['dendriteID'] + ' dendrite.',
    data = input['data'][sweeps],
    gain = 1.,
    unit = 'millivolt',
    electrode = input['electrode'],
    stimulus_description = 'N/A',
    timestamps = timestamps,
    comments = 'Somatic current clamp recording corresponding to ' +\
               'the initial part of the calcium imaging period ' +\
               'at ' + input['dendriteID'] + ' dendrite.' +\
               'The first dimension corresponds to individual ' +\
               'recording sweeps. The second dimension corresponds to ' +\
               'individual sweep data samples. The associated timestamps ' +\
               'variable provides timestamps for the second dimension. ' +\
               'This variable has to be combined with starting_time and ' +\
               'rate variables to get absolute timestamps ' +\
               'for each data point.   ' +\
               'data_continuity = step.   ' +\
               'rate = ' + str(input['imagingRate']))
  
  nwb.add_acquisition(current_clamp_series)
  return nwb


def reshapeData(data):
  '''
  reshapedArray = reshapeData(data)

  Function converts an ndarray of the shape (n, )
  into an ndarray of the shape (n, m, p).
  '''
  
  dim1size = len(data)
  dim2size = len(data[0])
  dim3size = len(data[0][0])
  reshapedArray = np.empty((dim1size,dim2size,dim3size))
  for element in range(dim1size):
    reshapedArray[element] = data[element]

  return reshapedArray