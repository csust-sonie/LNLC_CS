function [Y]=readbmp(bmpfile,begin,frames)
    
	path=['.\testbmp\' bmpfile '\'];
	fileFolder=fullfile(path);
	dirOutput=dir(fullfile(fileFolder,'*.bmp'));
	fileNames={dirOutput.name};
	sort(fileNames);

    for frame=1:frames
    	Y{frame}=double(rgb2gray(uint8(imread([path fileNames{frame+begin-1}]))));
    end

end