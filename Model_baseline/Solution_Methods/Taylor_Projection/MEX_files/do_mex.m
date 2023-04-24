folder_name='compiled_on_Windows10';

mkdir(folder_name);

for file = dir('*.F')'
    mex('-largeArrayDims', file.name, '-outdir',folder_name)
 end

