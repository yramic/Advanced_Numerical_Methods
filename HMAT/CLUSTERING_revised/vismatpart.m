function vismatpart(NF,FF)
% Visualization of near field and far field
    figure('Name','matrix partitioning');
    plot([0 1 1 0 0],[0 0 1 1 0],'k-'); hold on;
for nfbox = NF'
    plot([nfbox(1) nfbox(2)],[nfbox(3) nfbox(4)],'m*');
end
for ffbox = FF'
    plot([ffbox(1) ffbox(2) ffbox(2) ffbox(1) ffbox(1)],...
         [ffbox(3) ffbox(3) ffbox(4) ffbox(4) ffbox(3)],'b-');
end

end

