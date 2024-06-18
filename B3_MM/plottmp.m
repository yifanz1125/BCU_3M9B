close all   
for i=1:10
        figure(1);
        plot(Group.delta_stb(1:5e4,i),'Color',[(50+20*i)/255 150/255 (250-20*i)/255]); hold on;
%         strlb(no_group)="group"+no_group;
%         legend(strlb);
        cycle=size(Group.delta_stb(1:5e4,i),1);
        plot(Group.CUEP(i)*ones(cycle,1),'Color',[(50+20*i)/255 150/255 (250-20*i)/255],'LineStyle',':');    hold on;
%         strUEP(no_group)="CUEP"+no_group;
%         legend(strUEP);
%         figure(10+no_group);
%         plot(Group.omega_stb(:,no_group)); hold on;
%         figure(3)
%         plot(Group.delta_unstb(:,no_group)); hold on;
%         strlb(no_group)="group"+no_group;
%         legend(strlb);
%         plot(Group.CUEP(no_group)*ones(cycle_stb,1),':');    hold on;
%         strUEP(no_group)="CUEP"+no_group;
%         legend(strUEP);
    end