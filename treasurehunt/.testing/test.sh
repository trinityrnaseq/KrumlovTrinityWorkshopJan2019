#!/bin/bash

./encode_msg.pl --encode "LINE01 something cool
LINE02 whatever
LINE03 ultra cool" > ladeda.fa

./encode_msg.pl --decode ./ladeda.fa



