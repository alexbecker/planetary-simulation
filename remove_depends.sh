# Removes dependencies from the end of a Makefile

sed -i '/^.*\.o:.*$/d' Makefile
sed -i '/^ .*$/d' Makefile
