#!/bin/bash
# --------------------------------------------------------------------------- #
#   ==   ===  ==                                                              #
#  ||   ||=  || ))  support s. r. o. 2017, www.cfdsupport.com                 #
#   ==        ==                                                              #
# --------------------------------------------------------------------------- #
# Deleting Script: deletes fields "toDeleteFields" in processor* folders      #
#                  every "sleepHours" hours (for saving space in storage)     #
# --------------------------------------------------------------------------- #

toDeleteFields="U_0 Uar Uart Uat Uf Uf_0 Urt nuTilda_0 nuTilda nut pTot p_rgh phi"
sleepHours=2

while true; do
    echo Deleting fields $toDeleteFields
    for name in $toDeleteFields; do
        find ./processor* -name $name -not -wholename "*/0/*" | xargs rm -v
    done
    sleep ${sleepHours}h
done
