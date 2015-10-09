function s = getSNR()
    load lcmlog_2000_01_01_00.mat
    signal = attitude(:, 11);
    noise = -9.80665 - signal;
    s = snr(signal, noise);
end