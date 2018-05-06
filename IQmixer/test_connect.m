% dcSrcInterface =  ustc_dc_v1([49,50]);
% dcSrcObj = dcSource.GetInstance('myDCSource',dcSrcInterface);
% dcChnl1 = dcSrcObj.GetChnl(1);
% dcChnl2 = dcSrcObj.GetChnl(2);
% dcChnl1.dcval = 30000;
% dcChnl2.dcval = 30000;

% dcval1=-30000;
% dcval2=-30000;
% dcSrcInterface.SetDC(dcval1,1);
% dcSrcInterface.SetDC(dcval1,2);

%  daInterface = ustc_da_v1([49,50]);
%  awgObj = awg.GetInstance('myAWG',daInterface);
% daChnl1 = awgObj.GetChnl(1);
% daChnl2 = awgObj.GetChnl(2);