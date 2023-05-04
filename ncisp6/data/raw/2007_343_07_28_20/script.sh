old_prefix="NE"
old_suffix=".2007_343_07_28_20.BHR"

new_prefix="NE"
new_suffix='.r'


for file in $old_prefix*${old_suffix};do
	filename=$(basename "$file")
	echo $file
	prefix=${filename%${old_suffix}}
	new_filename="${prefix}${new_suffix}"
	cp "$file" "$new_filename"
done
