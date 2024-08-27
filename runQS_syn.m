clear;
nr = 20;
explist = 1:13;
nexp = length(explist);
qlist = 0.8:0.01:0.92;
p0 = 1;
option = 'Sto_block';
option_align = "L1";
for id = 1:nexp
    expid = explist(id);
    q = qlist(id);
    expstr = num2str(expid,'%03.f');
    testQS_syn
    save(strcat('syn_',option,num2str(expid),'.mat'),'p0','q','p1','p2','q1','q2'...
        ,'n','t_MPLS','t_LS','mean_error_LS_all','median_error_LS_all'...
        ,'t_LS5','mean_error_LS5_all','median_error_LS5_all'...
        ,'mean_error_MPLS_all','median_error_MPLS_all'...
        ,'mean_error_IRLS_all','median_error_IRLS_all','nb','idx');
    copyfile('testQS_syn.m',['testQS_syn_save' expstr '.m'],'f');
end