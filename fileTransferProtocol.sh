#!/bin/bash
# permanetly add IP-ADDRESS to ssh's known hosts for easy file transfer
# i.e. no password entry required
#ssh-keygen -t rsa
#ssh-copy-id -i rpi.pub USERNAME@IP-ADDRESS
# loop it to continuously send new temperature data written out 
# by the Amprobe TMD-56 thermometer
x=1
while [ $x -le 5 ]
do
    # find the newest file in the record folder
    # (run this .sh script in the record folder of TMD-56)
    # write endline to small file
    echo $(tail -1 $(ls -t | grep *REG*) | awk ' { print $4 } ') > ... 
    tempToTransfer
    # output new temperature in cmd
    cat tempToTransfer
    
    # transfer small file to RASPBERRY PI
    scp tempToTransfer USERNAME@IP-ADDRESS:/home/USERNAME/piHeater
    
    # arbitrary sleep time of 0.5 seconds
    sleep 0.5
done
