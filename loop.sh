find . -type f -name "*.dat" | while read -r dir; do
    # Remove each found file
    file="${dir#./}"
    echo ""
    echo $file
    git filter-repo --force --invert-paths --path-match $file
done
