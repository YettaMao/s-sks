old_prefix="NE"
old_suffix=".2007_288_12_29_33.BHT"

new_prefix="NE"
new_suffix='.t'


for file in $old_prefix*${old_suffix};do
	filename=$(basename "$file")
	echo $file
	prefix=${filename%${old_suffix}}
	new_filename="${prefix}${new_suffix}"
	cp "$file" "$new_filename"
done
