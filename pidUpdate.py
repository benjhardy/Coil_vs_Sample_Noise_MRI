# imports
from simple_pid import PID # open-source module
import numpy as np
import RPi.GPIO as GPIO
import time
import matplotlib.pyplot as plt
#import numpy as np
GPIO.setmode(GPIO.BOARD)
GPIO.setup(7,GPIO.OUT)
	
class controlled_system:
	"""Reads temp from the TMD series, turns on the relay if needed."""
	
def update(onOff,sv):
    
    # if 1 turn on relay, might want to return previous state
    if onOff < 0:
        GPIO.output(7,True)
        print('Turn on relay onOff, sv{}'.format((onOff,sv)))
    else:
        GPIO.output(7,False)
        print('Turn off relay onOff, sv{}'.format((onOff,sv)))
    #GPIO.cleanup()
    # reads temp from transferred data
    f = open('tempToTransfer',"r")
    v = float(f.readline())
    f.close()
    # read temp
    print('updated temp: {}'.format(v))
    return v
    
# Set PID parameters
pid = PID(1,.1,.05,setpoint=1)
pid.sample_time = 0.01 # in seconds update
pid.setpoint = 45 # celsius value to achieve
	
# Tune Kp,Ki,Kd
pid.tunings = (-1.35,0,6)

# declare class
v = controlled_system.update(23.0,pid.setpoint)

# init arrays for plotting
vals = np.array(v)
tempVals = np.array(v)

# update plot params
plt.ion()
t = 0
tv = np.array([t])

figure, ax = plt.subplots(figsize=(5,4))
line1, = ax.plot(t, vals,'--bo',label='PID val')
line2, = ax.plot(np.linspace(0,60,60),pid.setpoint*np.ones((60,)), ...
'--r', linewidth=2,label='Target Temp')
line3, = ax.plot(t,tempVals,'--ro',label='Thermocouple')

plt.title('PID tuning')
plt.xlabel('time points')
plt.ylabel('PID values')
plt.xlim([0,60])
plt.ylim([10,80])
plt.legend()

# update figure continuously to track PID and temperature
while True:
	# compute new output from the PID according to 
    #the systems current value
	onOff = pid(v)
	print('pid result: {}'.format(onOff))
 
	# feed the pid output to system and get the current value
	v = controlled_system.update(onOff,pid.setpoint)
	
	t += 1
	tv = np.append(tv,t)
	vals = np.append(vals,onOff+pid.setpoint)
	tempVals = np.append(tempVals,v)
	
	# updating data values
	line1.set_xdata(tv)
	line1.set_ydata(vals)
	line3.set_xdata(tv)
	line3.set_ydata(tempVals)

	# drawing updated values
	figure.canvas.draw()

	# This will run the GUI event
	# loop until all UI events
	# currently waiting have been processed
	figure.canvas.flush_events()

	time.sleep(0.01)
	
GPIO.cleanup()
