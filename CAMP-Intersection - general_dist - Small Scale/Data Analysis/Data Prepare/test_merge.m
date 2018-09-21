msn = [1:10,1:10]';
time = [1:20]';
sender_loc = [41:60]';
tt = table(time,msn,sender_loc);
rss = [101:120]';
received_idx = [1,2,5,8,12,17,18,19,20];

msn = msn(received_idx);
time = time(received_idx)+rand(length(received_idx),1);
rcd_loc = sender_loc(received_idx)+30;
rcd_rss = rss(received_idx);
rt = table(msn,time,rcd_loc,rcd_rss);
joined = outerjoin(tt,rt,'Keys','msn');

sortrows([joined(abs(joined.time_rt-joined.time_tt)<1,:);joined(isnan(joined.msn_rt),:)],'time_tt')