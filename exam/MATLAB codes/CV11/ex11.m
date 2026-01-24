clear;clc;close all;

%random solution
x = randi(2,1,20)-1
%evaluation
val = objective(x)

function val_sum = objective(x)
    val = zeros(45,1);
    weights = 1:45;
    val(1) = x(1) | not(x(2)) | not(x(20));
    val(2) = x(1) | x(17) | not(x(9));
    val(3) = not(x(1)) | x(2) | x(10);
    val(4) = x(2) | not(x(3)) | not(x(19));
    val(5) = x(3) | x(4) | not(x(5));
    val(6) = x(3) | not(x(7)) | x(11);
    val(7) = x(3) | x(9) | not(x(15));
    val(8) = x(4) | x(7) | not(x(16));
    val(9) = not(x(4)) | not(x(18)) | x(10);
    val(10) = x(4) | not(x(20)) | x(11);
    val(11) = x(5) | not(x(6));
    val(12) = not(x(5)) | not(x(14));
    val(13) = not(x(5)) | not(x(11));
    val(14) = not(x(5)) | x(6);
    val(15) = x(6) | x(7);
    val(16) = x(6) | x(14);
    val(17) = x(6) | not(x(15));
    val(18) = x(7) | x(20);
    val(19) = x(7) | not(x(11));
    val(20) = x(8) | x(13);
    val(21) = not(x(8)) | x(10);
    val(22) = x(9) | not(x(12));
    val(23) = not(x(9)) | not(x(10));
    val(24) = x(10) | x(15);
    val(25) = x(10) | x(16);
    val(26) = not(x(10)) | x(14);
    val(27) = x(11) | x(17);
    val(28) = not(x(11)) | x(20);
    val(29) = not(x(11)) | not(x(19));
    val(30) = x(12) | not(x(18));
    val(31) = x(12) | x(19);
    val(32) = x(12) | x(15);
    val(33) = not(x(13)) | not(x(14));
    val(34) = not(x(13)) | x(16);
    val(35) = x(13) | x(15);
    val(36) = x(14) | x(1);
    val(37) = x(15) | not(x(2));
    val(38) = x(16) | not(x(17)) | x(8);
    val(39) = x(17) | x(18) | not(x(19));
    val(40) = not(x(18)) | not(x(19)) | not(x(5));
    val(41) = not(x(10)) | x(11) | x(12);
    val(42) = not(x(6)) | x(20) | not(x(4));
    val(43) = not(x(2)) | x(19) | not(x(6));
    val(44) = not(x(14))| x(18) | not(x(2));
    val(45) = not(x(3)) | x(11) | x(14);
    val_sum = sum(weights*val);
end