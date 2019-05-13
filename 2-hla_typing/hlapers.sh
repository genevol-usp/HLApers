#!/bin/bash

usage() {
    echo "Usage: HLApers [options]"
    echo ""
    echo -e "\t-h | --help"
    echo -e "\t-m | --module"
    echo -e "\t-i | --index"
    echo -e "\t-fq1"
    echo -e "\t-fq2"
}

while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
    -h | --help)
	usage
	exit
	;;
    -m | --module)
	module="$2"
	shift 2
	;;
    -i | --index)
	index="$2"
	shift 2
	;;
    *)
	echo "ERROR: unknown parameter $1"
	usage
	exit 1
	;;
    esac
done




