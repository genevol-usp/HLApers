#/bin/bash

usage() {
    echo "Usage: HLApers [options]"
    echo ""
    echo "\t-h --help"
}

while [ "$1" != "" ]; do
    PARAM=`echo $1 | awk -F= '{print $1}'`
    VALUE=`echo $1 | awk -F= '{print $2}'`
    case $PARAM in
        -h | --help)
            usage
            exit
            ;;
        *)
            echo "ERROR: unknown parameter ${$PARAM}"
            usage
            exit 1
            ;;
    esac
    shift
done
