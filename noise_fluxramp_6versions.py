#v1,3,5,6 do not make big difference. For consistancy, use v3 for all later noise level comparison
#v1 and 3 are inverse in sign, 3 & 6 are the same. v1,3,6 are identical in PSD
#v5 is only a little different because it also includes the real part
#v4 is using the real part
#v2 is not an accountable method

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
	def __init__(self):
		self.I = []
		self.Q = []


	def readData(self, filename, theta0):
		datafile = open(filename)

		self.filetype = filename.split('.')[-1]

		newname = filename.split('/')[-1]
		self.ResName = filename.split('/')[-1].split('_')[0]
		self.att1 = filename.split('_')[-5]
		self.att2 = filename.split('_')[-3]

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

		if self.filetype =='iqfit':
			self.theta0 = -np.arctan2(self.Q_trans[self.FrInd], self.I_trans[self.FrInd]) #This is how the Matlab code rotates the circle
		else:
			self.theta0 = theta0

		self.I = (self.I_off - self.xc)*np.cos(self.theta0) - (self.Q_off - self.yc)*np.sin(self.theta0)
		self.Q = (self.I_off - self.xc)*np.sin(self.theta0) + (self.Q_off - self.yc)*np.cos(self.theta0)

		self.time_len = len(self.I)/self.fs
		self.time = np.arange(0, self.time_len, 1./self.fs)

		#Mag = np.sqrt(I**2 + Q**2)
		#Phase = np.arctan(Q/I)

		self.Mag = np.sqrt(self.I**2 + self.Q**2)
		self.Phase = np.arctan2(self.Q,self.I)

		sub = np.abs(self.F - self.Fr)
		self.FrInd = list(sub).index(min(sub))

		#plt.plot(self.Q[:300])
		#plt.plot(self.I[:300],self.Q[:300])
		#plt.show()
		#sys.exit()
		return self.theta0


	def plotIQ(self):

		Ind = self.FrInd
		#Ind = 98
		#if Current >= 0.0:
		#plt.plot(self.I[Ind], self.Q[Ind], linestyle = '', marker = '*', alpha = .5, label = 'att1='+att1+', att2='+att2)
		plt.plot(self.I[:200], self.Q[:200], linestyle = '-', marker = '', alpha = .5)#, label = 'Framp:%.0fKHz'%(f_ramp/1e3))
		#ax.plot(self.I_old[Ind], self.Q_old[Ind], linestyle = '', marker = '*', alpha = .5, markersize = 15)#, label = 'Framp:%.0fKHz'%(f_ramp/1e3))
		#ax.plot(self.I_old, self.Q_old, linestyle = '-', marker = 'o', alpha = .5)#, label = 'Framp:%.0fKHz'%(f_ramp/1e3))
		#ax.plot(self.I_off[Ind], self.Q_off[Ind], linestyle = '', marker = '*', alpha = .5)#, label = 'Framp:%.0fKHz'%(f_ramp/1e3))
		#ax.plot(self.I_off, self.Q_off, linestyle = '-', marker = '.', alpha = .5)#, label = 'Framp:%.0fKHz'%(f_ramp/1e3))

		plt.xlabel('Re [mV]', size = 30)
		plt.ylabel('Im [mV]', size = 30)

		#plt.savefig('i_vs_q_LO_%.1f_RF_%.1f.pdf'%(p1, p2))
		#plt.clf()


	def plotIQrel(self, len, color, islabel):
		Ind = self.FrInd
		if islabel=='on':
			plt.plot(self.I[:len]/self.r, self.Q[:len]/self.r, linestyle = '-', marker = '', alpha = 1, color = color, label = 'att1='+self.att1+', att2='+self.att2)
		else:
			plt.plot(self.I[:len]/self.r, self.Q[:len]/self.r, linestyle = '', marker = '.', alpha = .5, color = color)



		#plt.plot(self.I[Ind]/self.r, self.Q[Ind]/self.r, linestyle = '', marker = '*', alpha = .5)#, label = 'Framp:%.0fKHz'%(f_ramp/1e3))
		#ax.plot(self.I_old[Ind], self.Q_old[Ind], linestyle = '', marker = '*', alpha = .5, markersize = 15)#, label = 'Framp:%.0fKHz'%(f_ramp/1e3))
		#ax.plot(self.I_old, self.Q_old, linestyle = '-', marker = 'o', alpha = .5)#, label = 'Framp:%.0fKHz'%(f_ramp/1e3))
		#ax.plot(self.I_off[Ind], self.Q_off[Ind], linestyle = '', marker = '*', alpha = .5)#, label = 'Framp:%.0fKHz'%(f_ramp/1e3))
		#ax.plot(self.I_off, self.Q_off, linestyle = '-', marker = '.', alpha = .5)#, label = 'Framp:%.0fKHz'%(f_ramp/1e3))

		#plt.savefig('i_vs_q_LO_%.1f_RF_%.1f.pdf'%(p1, p2))
		#plt.clf()


	def demod(self, V_to_flux_slope):
		self.f_ramp = 10e3
		self.N_cycle = 2

		if V_to_flux_slope=='None':
			fc = self.N_cycle * self.f_ramp
		else:
			fc = 0

		N_ramp = int(self.time_len * self.f_ramp)
		N_pnt = self.fs/self.f_ramp
		self.time_new = np.arange(N_ramp-1)/self.f_ramp


		ind_begin = 0#int(0 * fs * 1e-6)
		ind_end = N_pnt#int(33 * fs * 1e-6)
		#print N_pnt
		#print 'ind_begin:ind_end', ind_begin, ind_end

		#plt.plot(self.Q[:300])
		#plt.show()
		#sys.exit()
		data = self.I +1j*self.Q
		data_initial = data[ind_begin:ind_end]
		#ratio = (data_initial[1:-1]/data_initial[0:-2]).mean()
		#fc_fit = np.arctan2(ratio.imag, ratio.real)/2/np.pi*fs



		Phase_v1 = np.empty(N_ramp-1)
		Phase_v2 = np.empty(N_ramp-1)
		Phase_v3 = np.empty(N_ramp-1)
		Phase_v4 = np.empty(N_ramp-1)
		Phase_v5 = np.empty(N_ramp-1)
		Phase_v6 = np.empty(N_ramp-1)

		self.Q_avg = 0.

		Mix = np.empty(N_ramp-1)
		for i in range(N_ramp-1):
			data_period = data[i*N_pnt+ind_begin : i*N_pnt+ind_end]
			time_period = self.time[i*N_pnt+ind_begin : i*N_pnt+ind_end]

			'''
			if i in [0,1,2, 1000,1001,1002, 1996,1997,1998]:
				plt.plot(data_period.imag, label=i)
				#plt.plot(data_period.real, data_period.imag)
			#plt.plot(data_period.imag)

			plt.plot(data_period.imag)
			plt.show()
			plt.exit()
			'''

			self.Q_avg += data_period.imag

			mix_v1 = (data_period.imag*np.exp(1j*2*np.pi*fc*time_period)).mean()
			phase_v1 = np.arctan2(mix_v1.imag, mix_v1.real)
			Phase_v1[i] = phase_v1

			mix_v2 = (data_period/data_initial).mean()
			phase_v2 = np.arctan2(mix_v2.imag, mix_v2.real)
			Phase_v2[i] = phase_v2

			#in v3, only the Q is used.
			mix_v3 = (data_period.imag/np.exp(1j*2*np.pi*fc*time_period)).mean()
			phase_v3 = np.arctan2(mix_v3.imag, mix_v3.real)
			Phase_v3[i] = phase_v3
			Mix[i] = mix_v3


			mix_v4 = (data_period.real/np.exp(1j*2*np.pi*fc*time_period)).mean()
			phase_v4 = np.arctan2(mix_v4.imag, mix_v4.real)
			Phase_v4[i] = phase_v4

			mix_v5 = (data_period/np.exp(1j*2*np.pi*fc*time_period)).mean()
			phase_v5 = np.arctan2(mix_v5.imag, mix_v5.real)
			Phase_v5[i] = phase_v5

			v6_fft = np.fft.fft(data_period.imag)
			mix_v6 = v6_fft[2]# - mix_fft[0]
			phase_v6 = np.arctan2(mix_v6.imag, mix_v6.real)
			Phase_v6[i] = phase_v6

		self.Q_avg /= N_ramp-1

		self.Phase_v1 = Phase_v1 - Phase_v1[0]
		self.Phase_v2 = Phase_v2 - Phase_v2[0]
		self.Phase_v3 = Phase_v3 - Phase_v3[0]
		self.Phase_v4 = Phase_v4 - Phase_v4[0]
		self.Phase_v5 = Phase_v5 - Phase_v5[0]
		self.Phase_v6 = Phase_v6 - Phase_v6[0]

		#self.Phase_v1 = signal.detrend(Phase_v1, type = 'linear')#'constant')
		#self.Phase_v2 = signal.detrend(Phase_v2, type = 'linear')#'constant')
		self.Phase_v3 = signal.detrend(Phase_v3, type = 'linear')#'constant')
		#self.Phase_v4 = signal.detrend(Phase_v4, type = 'linear')#'constant')
		#self.Phase_v5 = signal.detrend(Phase_v5, type = 'linear')#'constant')
		#self.Phase_v6 = signal.detrend(Phase_v6, type = 'linear')#'constant')

		if V_to_flux_slope=='None':
			self.Phase = self.Phase_v3
		else:
			self.Phase = Mix-Mix[0]



	def PhaseToV(self):
		
		def funcWave(x, Vamp, a, b, v0, ratio): #ratio = Ls/Lj
			m = a*x + b #phase
			y = Vamp*np.sin(m)/(1+ratio*np.sin(m)) + v0
			#y = Vamp*np.cos(m)/(1+ratio*np.cos(m)) + v0
			return y

		self.phase_ramp = np.arange(0, self.N_cycle * 2*np.pi, self.N_cycle * 2*np.pi/len(self.Q_avg))
		p0 = np.array([0.1, 1., 0.1, 0, 0.4])
		params = opt.curve_fit(funcWave, self.phase_ramp, self.Q_avg, p0)
		self.Vamp = params[0][0]
		self.a = params[0][1]
		self.b = params[0][2]
		self.v0 = params[0][3]
		self.ratio = params[0][4]

		print 'Vamp=',self.Vamp
		print 'a=',self.a
		print 'b=',self.b
		print 'v0=',self.v0
		print 'ratio=',self.ratio


		self.phase_fit = np.arange(0, self.N_cycle * 2*np.pi, 1e-3)
		self.V_fit = funcWave(self.phase_fit, self.Vamp, self.a, self.b, self.v0, self.ratio)

		SlopeInd = list(abs(self.Q_avg)).index(min(abs(self.Q_avg-0)))

		self.V_to_flux = np.diff(self.Q_avg)/np.diff(self.phase_ramp/(2*np.pi))
		self.V_to_flux_slope = self.V_to_flux[SlopeInd]

		print 'V_to_flux_slope: ', self.V_to_flux_slope

		return self.V_to_flux_slope

	def plotDemodIQ(self):
		plt.plot((self.a*self.phase_fit+self.b)/(2*np.pi), self.V_fit, color='r')
		plt.plot((self.a*self.phase_ramp+self.b)/(2*np.pi), self.Q_avg, linestyle='', marker='o', markerfacecolor='None', markeredgecolor='b')
		plt.plot(self.phase_ramp/(2*np.pi), self.Q_avg, linestyle='', marker='.', color='g', alpha = 0.5)
		plt.plot(-self.Phase/(2*np.pi), np.array([0]*len(self.Phase)), linestyle='', marker='.')




	def plotDemod(self):
		plt.plot(self.time_new, self.Phase_v1, color = 'r', label = 'v1', alpha = .3)
		plt.plot(self.time_new, self.Phase_v2, color = 'y', label = 'v2', alpha = .3)
		plt.plot(self.time_new, self.Phase_v3, color = 'b', label = 'v3', alpha = .3)
		plt.plot(self.time_new, self.Phase_v4, color = 'g', label = 'v4', alpha = .3)
		plt.plot(self.time_new, self.Phase_v5, color = 'm', label = 'v5', alpha = .3)
		plt.plot(self.time_new, self.Phase_v6, color = 'c', label = 'v6', alpha = .3)
		plt.ylabel('phase [rad]')
		plt.legend(loc = 'upper left')
		plt.show()



	def mkPhase(self):
		f = open(nsfilename.strip('ns')+'phase', 'w')
		for i in range(len(Phase_v3)):
			f.write(str(Phase_v3[i]) + '\n')
		f.close
		sys.exit()

	def calPSD(self, V_to_flux_slope):
		#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		#calculate psd
		#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		#window = 'hanning'
		window = 'boxcar' #based on the old Matlab code

		if V_to_flux_slope=='None':
			fs = self.f_ramp
			nperseg = 500

			self.flux = self.Phase/2./np.pi
			self.ff_ns, self.flux_psd = signal.welch(self.flux, fs = fs, window = window, nperseg = nperseg, noverlap = 0, return_onesided=True)
		
		else:
			fs = 1e4
			nperseg = 500
			self.flux = self.Phase/ V_to_flux_slope
			self.ff_ns, self.flux_psd = signal.welch(self.flux, fs = fs, window = window, nperseg = nperseg, noverlap = 0, return_onesided=True)


		'''
		flux_v1 = self.Phase_v1/2./np.pi
		self.ff_ns_v1, self.flux_psd_v1 = signal.welch(flux_v1, fs = fs, window = window, nperseg = nperseg, noverlap = 0, return_onesided=True)

		flux_v2 = self.Phase_v2/2./np.pi
		self.ff_ns_v2, self.flux_psd_v2 = signal.welch(flux_v2, fs = fs, window = window, nperseg = nperseg, noverlap = 0, return_onesided=True)

		flux_v3 = self.Phase_v3/2./np.pi
		self.ff_ns_v3, self.flux_psd_v3 = signal.welch(flux_v3, fs = fs, window = window, nperseg = nperseg, noverlap = 0, return_onesided=True)

		flux_v4 = self.Phase_v4/2./np.pi
		self.ff_ns_v4, self.flux_psd_v4 = signal.welch(flux_v4, fs = fs, window = window, nperseg = nperseg, noverlap = 0, return_onesided=True)

		flux_v5 = self.Phase_v5/2./np.pi
		self.ff_ns_v5, self.flux_psd_v5 = signal.welch(flux_v5, fs = fs, window = window, nperseg = nperseg, noverlap = 0, return_onesided=True)

		flux_v6 = self.Phase_v6/2./np.pi
		self.ff_ns_v6, self.flux_psd_v6 = signal.welch(flux_v6, fs = fs, window = window, nperseg = nperseg, noverlap = 0, return_onesided=True)
		'''

		#Phase_v6_detrend = signal.detrend(Phase_v6, type = 'linear')
		#flux_v6_detrend = Phase_v6_detrend/2./np.pi
		#ff_ns_v6_detrend, flux_psd_v6_detrend = signal.welch(flux_v6_detrend, fs = fs, window = window, nperseg = nperseg, noverlap = 0, return_onesided=True)


	def plotPSD(self, V_to_flux_slope, phaseAxis, psdAxis):
		#fig = plt.figure()

		#phaseAxis = fig.add_subplot(211)
		if V_to_flux_slope=='None':
			phaseAxis.plot(self.time_new, self.flux*1e3, label = 'flux-ramp demod', alpha=.5, color='b')#'att1='+self.att1+', att2='+self.att2)
		else:
			phaseAxis.plot(self.time_new, self.flux*1e3, label = 'open-loop', alpha=.5, color='g')#'att1='+self.att1+', att2='+self.att2)


		#phaseAxis.plot(self.time_new, self.Phase_v1, label = 'v1')
		#phaseAxis.plot(self.time_new, self.Phase_v2, label = 'v2')
		#phaseAxis.plot(self.time_new, self.Phase_v3, label = 'v3')
		#phaseAxis.plot(self.time_new, self.Phase_v4, label = 'v4')
		#phaseAxis.plot(self.time_new, self.Phase_v5, label = 'v5')
		#phaseAxis.plot(self.time_new, self.Phase_v6, label = 'v6')
		#phaseAxis.plot(self.time_new, Phase_v6_detrend, label = 'from v6 detrend')
		#phaseAxis.plot(self.time_new, flux_v6_detrend, label = 'from v6 detrend')
		phaseAxis.set_xlabel('Time [s]')
		#phaseAxis.set_ylabel('Phase [rad]')
		phaseAxis.set_ylabel('m$\Phi_0$')
		phaseAxis.legend(loc = 'lower left')
		#plt.subplots_adjust(left=0.1, right=0.7, top=0.9, bottom=0.1)
		#plt.legend(bbox_to_anchor=(1.5, 1.))

		#psdAxis = fig.add_subplot(212)
		if V_to_flux_slope=='None':
			psdAxis.plot(self.ff_ns, np.sqrt(self.flux_psd)*1e6, alpha=.5, color='b')
		else:
			psdAxis.plot(self.ff_ns, np.sqrt(self.flux_psd)*1e6, alpha=.5, color='g')

		#psdAxis.plot(self.ff_ns_v1, np.sqrt(self.flux_psd_v1)*1e6, label = 'v1')
		#psdAxis.plot(self.ff_ns_v2, np.sqrt(self.flux_psd_v2)*1e6, label = 'v2')
		#psdAxis.plot(self.ff_ns_v3, np.sqrt(self.flux_psd_v3)*1e6, label = 'v3')
		#psdAxis.plot(self.ff_ns_v4, np.sqrt(self.flux_psd_v4)*1e6, label = 'v4')
		#psdAxis.plot(self.ff_ns_v5, np.sqrt(self.flux_psd_v5)*1e6, label = 'v5')
		#psdAxis.plot(self.ff_ns_v6, np.sqrt(self.flux_psd_v6)*1e6, label = 'v6')
		psdAxis.set_xlabel('Freq [Hz]')
		psdAxis.set_ylabel('Flux Noise [$\mu \Phi_0/\sqrt{Hz}$]')
		psdAxis.set_xscale('log')
		psdAxis.set_yscale('log')
		psdAxis.set_ylim(0.01, 1e4)
		#psdAxis.legend(loc = 'lower left')


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#read in noise data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#noifilelist = glob.glob('../../Muxchip/03-17-2016_IQBox/*.ns')
#iqfilelist = glob.glob('../../Muxchip/03-17-2016_IQBox/*.iqfit')
#noifilelist = glob.glob('../10-10-2016_IQBox/flux-ramp_noise/*.ns')
#iqfilelist = glob.glob('../10-13-2016_IQBox/flux-ramp_noise_55uA/*.iqfit')
#flnoifilelist = glob.glob('../10-13-2016_IQBox/flux-ramp_noise_55uA/*.ns')
opnoifilelist = glob.glob('../10-12-2016_IQBox/open-loop_noise_55uA/*.ns')
flnoifilelist = glob.glob('../10-12-2016_IQBox/flux-ramp_noise_55uA/*.ns')
iqfilelist = glob.glob('../10-12-2016_IQBox/flux-ramp_noise_55uA/*.iqfit')

cmap = ['r','y','g','b','m','c','k']

for i in [1]:#[0,3,6,9,12,15]:
	iqfilename = iqfilelist[i]
	flnoifilename = flnoifilelist[i]
	opnoifilename = opnoifilelist[i]
	print iqfilename
	print flnoifilename
	print opnoifilename

	color = cmap[i%7]

	ResIQ = Res()
	theta0= ResIQ.readData(iqfilename, 0)
	#ResIQ.plotIQ()
	#ResIQ.plotIQrel(200, color, 'on')

	fig = plt.figure()
	phaseAxis = fig.add_subplot(211)
	psdAxis = fig.add_subplot(212)

	
	flNoiIQ = Res()
	flNoiIQ.readData(flnoifilename, theta0)
	#flNoiIQ.plotIQ()
	#flNoiIQ.plotIQrel(200, color, 'off')
	flNoiIQ.demod('None')
	V_to_flux_slope = flNoiIQ.PhaseToV()
	#NoiIQ.plotDemodIQ()
	#NoiIQ.plotDemod()
	flNoiIQ.calPSD('None')
	flNoiIQ.plotPSD('None', phaseAxis, psdAxis)
	

	opNoiIQ = Res()
	opNoiIQ.readData(opnoifilename, theta0)
	#opNoiIQ.plotIQ()
	opNoiIQ.demod(V_to_flux_slope)
	opNoiIQ.calPSD(V_to_flux_slope)
	opNoiIQ.plotPSD(V_to_flux_slope, phaseAxis, psdAxis)
	plt.ylim(1, 1e2)

#plt.legend(loc='upper right')
plt.tight_layout()
plt.suptitle(opNoiIQ.ResName)
plt.show()
