
import numpy as np
from scipy import signal
import struct;
import sys, os
import glob
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rc('figure', facecolor = 'w')


class Res():
	#define a class for data handling
	def __init__(self):
		self.I = []
		self.Q = []


	def readIQ(self, filename):
		#read in the .iqfit file, return trans-rotate parameters
		self.filename=filename
		datafile = open(filename)

		self.filetype = filename.split('.')[-1]

		newname = filename.split('/')[-1]
		self.Temp = 100
		#self.Temp = float(newname.split('_')[5].strip('K'))*1e3
		self.ResName = filename.split('/')[-1].split('_')[0]
		self.att1 = newname.split('_')[2]
		self.att2 = newname.split('_')[4]

		data = []
		for line in datafile:
		    data.append(line)
		datafile.close()

		self.Fr = float(data[5].split()[1])
		self.Fr *= 1e9 #fitted resonance frequency

		self.x_offset = float(data[3].split()[1])
		self.y_offset = float(data[4].split()[1])

		self.xc = float(data[6].split()[1])
		self.yc = float(data[7].split()[1])
		self.r = float(data[8].split()[1])
		self.theta0 = float(data[12].split()[1])

		self.fs = float(data[27].split()[1])


		x = []
		y = []
		f = []
		for i in range(35, len(data)):
			f.append(float(data[i].split()[0]))
			x.append(float(data[i].split()[1]))
			y.append(float(data[i].split()[2]))
		self.F = np.array(f)*1e9
		self.I_old = np.array(x)
		self.Q_old= np.array(y)

		sub = np.abs(self.F - self.Fr)
		self.FrInd = list(sub).index(min(sub))

		delay = 30e-9 #cable delay
		self.phase_delay = -self.F*delay*2*np.pi

		
		self.I_off = (self.I_old - self.x_offset)*np.cos(self.phase_delay) - (self.Q_old - self.y_offset)*np.sin(self.phase_delay)
		self.Q_off = (self.I_old - self.x_offset)*np.sin(self.phase_delay) + (self.Q_old - self.y_offset)*np.cos(self.phase_delay)

		self.I_trans = self.I_off - self.xc
		self.Q_trans = self.Q_off - self.yc

		self.theta0 = -np.arctan2(self.Q_trans[self.FrInd], self.I_trans[self.FrInd]) #This is how the Matlab code rotates the circle

		self.I = (self.I_off - self.xc)*np.cos(self.theta0) - (self.Q_off - self.yc)*np.sin(self.theta0)
		self.Q = (self.I_off - self.xc)*np.sin(self.theta0) + (self.Q_off - self.yc)*np.cos(self.theta0)

		self.time_len = len(self.I)/self.fs
		self.time = np.arange(0, self.time_len, 1./self.fs)

		self.Mag = np.sqrt(self.I**2 + self.Q**2)
		self.Phase = np.arctan2(self.Q,self.I)

		sub = np.abs(self.F - self.Fr)
		self.FrInd = list(sub).index(min(sub))


		return self.x_offset, self.y_offset, -self.Fr*delay*2*np.pi, self.xc, self.yc, self.theta0

	def readNoi(self, filename, x_offset, y_offset, phase_delay, xc, yc, theta0):
		#read in the noise file or pulse file, do some trans-rotation based on the .iqfit file
		#after trans-rotation, the Q part is used for demodulation
		self.filename=filename
		datafile = open(filename)
		newname = filename.split('/')[-1]

		x = []
		y = []
		for line in datafile:
			try:
				x.append(float(line.split()[0]))
				y.append(float(line.split()[1]))
			except ValueError:
				continue
		self.I_old = np.array(x)
		self.Q_old = np.array(y)


		self.I_off = (self.I_old - x_offset)*np.cos(phase_delay) - (self.Q_old - y_offset)*np.sin(phase_delay)
		self.Q_off = (self.I_old - x_offset)*np.sin(phase_delay) + (self.Q_old - y_offset)*np.cos(phase_delay)

		self.I_trans = self.I_off - xc
		self.Q_trans = self.Q_off - yc


		self.I = (self.I_off - xc)*np.cos(theta0) - (self.Q_off - yc)*np.sin(theta0)
		self.Q = (self.I_off - xc)*np.sin(theta0) + (self.Q_off - yc)*np.cos(theta0)


		self.fs = 500e3 #sampling frequency of the digitizer
		self.time_len = len(self.I)/self.fs #how long the data stream is
		self.time = np.arange(0, self.time_len, 1./self.fs) #generate a time stamp


	def plotIQ(self):
		#a plottig function
		plt.plot(self.I[:200], self.Q[:200], linestyle = '-', marker = '', alpha = .5)#, label = 'Framp:%.0fKHz'%(f_ramp/1e3))

		plt.xlabel('Re [mV]', size = 30)
		plt.ylabel('Im [mV]', size = 30)

		#plt.savefig('i_vs_q_LO_%.1f_RF_%.1f.pdf'%(p1, p2))
		#plt.clf()

	def plotIQrel(self, len, color, islabel):
		#a plotting function
		Ind = self.FrInd
		if islabel=='on':
			plt.plot(self.I[:len]/self.r, self.Q[:len]/self.r, linestyle = '-', marker = '', alpha = 1, color = color, label = 'att1='+self.att1+', att2='+self.att2)
		else:
			plt.plot(self.I[:len]/self.r, self.Q[:len]/self.r, linestyle = '', marker = '.', alpha = .5, color = color)


	def demod(self):
		#the demodulation function, you may hate how I name the variables
		self.f_ramp = 10e3 #ramp frequency
		self.N_cycle = 2 #number of cycles in one ramp
		fc = self.N_cycle * self.f_ramp #carrier frequency

		N_ramp = int(self.time_len * self.f_ramp) #number of ramps in this data stream
		N_pnt = self.fs/self.f_ramp #number of points in one ramp, the demode window size will be determined based on this
		self.time_new = np.arange(N_ramp)/self.f_ramp #generate a new time stamp for the demode data
		ind_begin = 0 #window beginning index
		ind_end = N_pnt #window end index, here I am not throwing data points
		Phase_v3 = np.empty(N_ramp) #initialize an empty array to save the demode phase

		for i in range(N_ramp):
			#demode sectiong by section, a hamming window is applied to the ramp section
			Q_period = self.Q[i*N_pnt+ind_begin : i*N_pnt+ind_end]*np.hamming(ind_end - ind_begin)
			time_period = self.time[i*N_pnt+ind_begin : i*N_pnt+ind_end]

			#the exact demodulation formula, I name it "mix", because it is like the software version of an IQ mixer. "v3" is kept for documentary purpose. I tried 6 version.
			mix_v3 = (Q_period/np.exp(1j*2*np.pi*fc*time_period)).mean() 

			#out of the demodulated data, calculate the phase by arctan
			phase_v3 = np.arctan2(mix_v3.imag, mix_v3.real)

			#append the demode phase to the list
			Phase_v3[i] = phase_v3

		self.Phase = Phase_v3

	def twoPi(self):
		#this function removes the 2pi jump
		#it may only work for this specific type of data. It's a non-intelligent & violent function
		diff = np.diff(self.Phase)
		a = []
		b = []
		for i in range(len(diff)):
			if diff[i] > 0.2*2*np.pi:
				a.append(i+1)
			elif diff[i] < -0.8*2*np.pi:
				b.append(i+1)
			else:
				pass

		for j in range(len(a)):
			self.Phase[a[j]:b[j]] -= 2*np.pi

		#2pi in phase corresponds to one quantum flux
		self.flux = self.Phase/2./np.pi


	def plotDemod(self):
		#a plotting function
		plt.plot(self.time_new, self.Phase, color = 'b', label = 'v3', alpha = .3)
		plt.ylabel('phase [rad]')
		plt.legend(loc = 'upper left')
		plt.show()

	def writeFlux(self):
		#write the demode data in file	
		f = open(self.filename.strip('ns')+'flux', 'w')
		for i in range(len(self.flux)):
			f.write(str(self.time_new[i])+'\t'+str(self.flux[i]) + '\n')
		f.close()
	
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#main function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def Exct(iqfilename, dirname, res):
	flnoifilelist = glob.glob(dirname+res+'*.ns')
	print iqfilename
	Num = len(flnoifilelist)
	for i in range(7, 8):
		flnoifilename = flnoifilelist[i]
		print flnoifilename

		#read in the .iqfit file and return the transrotation parameters
		ResIQ = Res()
		x_offset, y_offset, phase_delay, xc, yc, theta0 = ResIQ.readIQ(iqfilename)
		
		#read in the noise or pulse file, and demodulate
		flNoiIQ = Res()
		flNoiIQ.readNoi(flnoifilename, x_offset, y_offset, phase_delay, xc, yc, theta0)
		flNoiIQ.demod()
		flNoiIQ.twoPi()

		#flNoiIQ.writeFlux()




resList = ['Res-8']
for res in resList:
	dirname = '../12-16-2016/pulse_files/'
	iqfilename = glob.glob(dirname+res+'_Att1_10_Att2_40_0.075K_Circle.iqfit')[0]
	Exct(iqfilename, dirname, res)
