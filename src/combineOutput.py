#!/usr/bin/python
# PYTHON script to create a global XDMF file from the local XDMF components

print ('---------------------------------------------------')
print (' Combining local XDMF files to produce global grid ')
print ('---------------------------------------------------')

import os

# create the header portion of the global XDMF file

output = open('../out/global_header.xmf', 'w')
output.write('<?xml version=\"1.0\" ?>\n')
output.write('<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n')
output.write('<Xdmf version = \"2.0\">\n')
output.write('<Domain>\n')
output.write('<Grid Name="time series" GridType="Collection" CollectionType="Temporal">\n')
output.close()

os.system('cat ../out/global_header.xmf > ../out/time_series.xmf')

num_frame  = 10
frame_rate = 100

for time in range (0,num_frame+1):
  frame = time*frame_rate
  print ('\nAssembling XDMF output for frame = %d\n' % frame)

  output = open('../out/loopxmf.xmf','w')
  output.write('<Grid Name="global mesh" GridType="Collection" CollectionType="Spatial">\n')
  output.write('<Time Value = "%d" />\n' % frame)
  output.close()

  os.system('cat ../out/loopxmf.xmf >> ../out/time_series.xmf')
  # cocatenate XDMF files from individual processors to create data for this time level
  filename = 'ls -l ../out/xdmf/data_' + str(frame).zfill(6) + '*'
  os.system(filename)
  filename = 'cat ../out/xdmf/data_' + str(frame).zfill(6) + '* >> ../out/time_series.xmf'
  os.system(filename)
  output = open('../out/loopxmf.xmf','w')
  output.write('</Grid>\n')
  output.close()
  os.system('cat ../out/loopxmf.xmf >> ../out/time_series.xmf')

  # create the footer portion of the global XDMF file for this time level

output = open('../out/global_footer.xmf', 'w')
output.write('</Grid> <!-- time series ends -->\n')
output.write('</Domain>\n')
output.write('</Xdmf>\n')
output.close()

# combine the header, main section and footer to create the global XDMF file

os.system('cat ../out/global_footer.xmf >> ../out/time_series.xmf')

# remove the header, main and footer files because we do not need them anymore

os.system('rm ../out/global_header.xmf')
os.system('rm ../out/loopxmf.xmf')
os.system('rm ../out/global_footer.xmf')

print ('Done')
