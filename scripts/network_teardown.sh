sudo tc qdisc del dev $DEV root

unset DEV
unset SPEED
unset BURST
unset LAT