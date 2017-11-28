import numpy as np

def readLithoDat(datFpath):

	with open(datFpath,"r") as ins:
		lines=[]
		for line in ins:
			lines.append(line)

	x=[];y=[];z=[];litho=[] 
	for i in range(len(lines)):
		if len(lines[i])>2: # last few lines may be empty
			tmp=lines[i].split()
			x.append(tmp[0]);y.append(tmp[1])
			z.append(tmp[2]);litho.append(tmp[3])
	return x,y,z,litho

def writeWellSgemsFileLitho(well,x,y,z,litho,slowness):

	fname=well+'.dat';fid=open(fname,'w')
	fid.write('Well '+well[-3:]+'\n'+'4'+'\n'+'x'+'\n'+'y'+'\n'+'z'+'\n'+'lithology'+'\n')
	for i in range(len(z)):
		fid.write(x[i]+' '+y[i]+' '+z[i]+' '+litho[i]+'\n')
	fid.close()
	
def writeWellSgemsFileSlow(well,x,y,z,litho,slowness):

	fname=well+'slow-ABV.dat';fid=open(fname,'w')
	fid.write('Well '+well[-3:]+'\n'+'4'+'\n'+'x'+'\n'+'y'+'\n'+'z'+'\n'+'slowness'+'\n')
	for i in range(350):#len(z)):
		fid.write(x[i]+' '+y[i]+' '+z[i]+' '+str(slowness[-int(litho[i])])+'\n')
	fid.close()

slowMin,slowMax=1./3228,1./1511 # from queen et. al
numLitho=5#12
slowness=np.linspace(slowMin,slowMax,numLitho)	

writeAllFilesFull=False
if writeAllFilesFull:
	well='Well18-'
	datFpath='../18_1_Lithos.dat'
	x,y,z,litho=readLithoDat(datFpath);
	writeWellSgemsFileSlow(well,x,y,z,litho,slowness)

	well='Well18A'
	datFpath='../18A_1_Litho.dat'
	x,y,z,litho=readLithoDat(datFpath)
	writeWellSgemsFileSlow(well,x,y,z,litho,slowness)

	well='Well18B'
	datFpath='../18B_31_Litho.dat'
	x,y,z,litho=readLithoDat(datFpath)
	writeWellSgemsFileSlow(well,x,y,z,litho,slowness)

	well='Well27-'
	datFpath='../27_1_Litho.dat'
	x,y,z,litho=readLithoDat(datFpath)
	writeWellSgemsFileSlow(well,x,y,z,litho,slowness)

	well='Well46-'
	datFpath='../46_1_Litho.dat'
	x,y,z,litho=readLithoDat(datFpath)
	writeWellSgemsFileSlow(well,x,y,z,litho,slowness)

	well='Well46A'
	datFpath='../46A_1_Litho.dat'
	x,y,z,litho=readLithoDat(datFpath)
	writeWellSgemsFileSlow(well,x,y,z,litho,slowness)

	well='Well47A'
	datFpath='../47A_1_Litho.dat'
	x,y,z,litho=readLithoDat(datFpath)
	writeWellSgemsFileSlow(well,x,y,z,litho,slowness)

	well='Well47C'
	datFpath='../47C_1_Litho.dat'
	x,y,z,litho=readLithoDat(datFpath)
	writeWellSgemsFileSlow(well,x,y,z,litho,slowness)

	well='Well55-'
	datFpath='../55_1_Litho.dat'
	x,y,z,litho=readLithoDat(datFpath)
	writeWellSgemsFileSlow(well,x,y,z,litho,slowness)

	well='Well56-'
	datFpath='../56_1_Litho.dat'
	x,y,z,litho=readLithoDat(datFpath)
	writeWellSgemsFileSlow(well,x,y,z,litho,slowness)

	well='Well56A'
	datFpath='../56A_1_Litho.dat'
	x,y,z,litho=readLithoDat(datFpath)
	writeWellSgemsFileSlow(well,x,y,z,litho,slowness)

	well='Well56B'
	datFpath='../56B_1_Litho.dat'
	x,y,z,litho=readLithoDat(datFpath)
	writeWellSgemsFileSlow(well,x,y,z,litho,slowness)

	well='Well57-'
	datFpath='../57_1_Litho.dat'
	x,y,z,litho=readLithoDat(datFpath)
	writeWellSgemsFileSlow(well,x,y,z,litho,slowness)
		
	well='Well15-'
	datFpath='../15_12RD_litho.dat'
	x,y,z,litho=readLithoDat(datFpath)
	writeWellSgemsFileSlow(well,x,y,z,litho,slowness)

	
well='Well56A'
datFpath='../56A_1_Litho.dat'
x,y,z,litho=readLithoDat(datFpath)
writeWellSgemsFileSlow(well,x,y,z,litho,slowness)



