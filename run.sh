find . -executable -type f | while read line; do
    ./$line
done
