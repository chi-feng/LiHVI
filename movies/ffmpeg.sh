#!/bin/bash
export prefix=$(zenity --entry --text "Prefix")

echo "Convering ppm to jpg"
for pic in *.ppm
do
pnmtojpeg "${pic}" > "${pic/%ppm/jpg}"
done

ffmpeg -r 60 -b 10000k -i ${prefix}.%04d.jpg ${prefix}1800.mp4
if zenity --question --text="Delete jpeg files?"; then
	rm ${prefix}*.jpg
	echo "Done."
else
	echo "Keeping jpeg files"
fi
if zenity --question --text="Delete ppm files?"; then
        rm ${prefix}*.ppm
        echo "Done."
else
        echo "Keeping ppm files"
fi



