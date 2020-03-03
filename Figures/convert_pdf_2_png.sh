

for i in *.pdf; do sips -s format png ${i} --out ${i/.pdf/.png}; done

