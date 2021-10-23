function plt_coh_temp(titl,GestList,nf,coh_temp)   % args: connect_matrix, component
    figure;
    plot(nf,coh_temp);
    hold on
    plot(nf,mean(coh_temp,2),'k','linew',1.5);
    plot(nf,median(coh_temp,2),'r','linew',1.5);
    set(gca,'xlim',[0 40]);
    legend([GestList,'Mean','Median']);
    title(strcat("MSC - ",titl));
end