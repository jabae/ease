function reader = getReader(key)

filename = fetch1(nda.Scan & key, 'filename');
path = 'T:\Two-Photon\Jake\160308';
fprintf('Loading from %s\n', path);
path = getLocalPath(fullfile(path, sprintf('%s_*.tif', filename)));

reader = ne7.scanimage.Reader5(path);

end