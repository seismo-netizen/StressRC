clear; clc; close all;
%% merge figures generated by step 1 to avi video
% Haoran Meng, Jan 22, 2020

% Mainshock Mw 7.1
writerObj = VideoWriter('./RC71_PK_STRESS.avi'); flag = 1;
% Foreshock Mw 6.4
% writerObj = VideoWriter('./RC64_PK_STRESS.avi'); flag = 2;
writerObj.FrameRate = 20;
open(writerObj);

for k = 1:200
    
    if flag == 1
        if k < 10
            fnm = ['./figures/RC71_PK_STRESS_00' num2str(k) '.png' ];
        elseif k<100 && k>=10
            fnm = ['./figures/RC71_PK_STRESS_0' num2str(k) '.png'];
        else
            fnm = ['./figures/RC71_PK_STRESS_' num2str(k) '.png' ];
        end
    end
    
    if flag == 2
        if k < 10
            fnm = ['./figures/RC64_PK_STRESS_00' num2str(k) '.png' ];
        elseif k<100 && k>=10
            fnm = ['./figures/RC64_PK_STRESS_0' num2str(k) '.png' ];
        else
            fnm = ['./figures/RC64_PK_STRESS_' num2str(k) '.png' ];
        end
    end

    thisimage = imread(fnm);
    writeVideo(writerObj, thisimage);
end
close(writerObj);
