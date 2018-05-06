% author:F.Liang
% data:2018/4/1
% version:1.0
% filename:FTDAs.m
% describe:�����ﶨ�ƿ�߾���DCԴ�����࣬��װ�棬ÿ�ζ�д���п���tcp���ӡ�
% ĿǰDA������tcpip serverģʽ�����ܵ�ָ����ʽֻ��һ�֣��ַ�������'DA=1;RW=1;ADDR=0x01;VAL=0x00000;'
% DA=(1,2,3,4);RW=(1W,0R);ADDR=0x(01,02,03,04);VAL=0x(00000~FFFFF);
% DA��4��ͨ�����豸��ÿ��ͨ����Ӧһ�Ի������������ԭ��������������ȶ��Ի�ϸ�������ã���Ϊ�������ⲿ����������Ԫ��������Ӱ�����ܡ�
% DA�ڽ���RW=1;��д����ʱ��û���κη���ֵ���ڽ���RW=0;ʱ��ֱ�ӷ��ش�DACָ��ͨ���������ĵ�ѹ����ֵ������������ʽΪ�ַ���'0x00000\n'?�������ɺ��ڵ���Ӳ���޸�����������ơ�
% RW=0ʱ��VAL�����壬�����ǰ��淶��ʽ����ȥ�����Դ���0���롣
% ADDR�����û�ֻʹ��01��������ַΪDAC�������������ͼĴ�����ֻ�Ը߼��û�ʹ�ã����Գ����д��������Ҫ�û��ṩ��ַ��Ϣ��������д����
% VALΪ��ѹ����ֵ��20bit�������룬0x00000��С����Ӧ�����Χ�п�����-2V~+2V��-3.5V~+3.5V��-7V~+7V���ӵ�·ֱ���������������������ڣ���Ҫ���ܷ�ǣ���Ŀǰ�汾Ϊ-2V~+2V
% �����豸���ܲ���18bitDAC������ֵ��ʽ���䣬��20bit��DAC��ֻ�������bit�����塣Ŀǰ�汾Ϊ20bit��
% ʹ�÷�����
% dac_ip='10.0.200.1';
% dac = FTDAs(dac_ip);
% dac.SetValue(1,hex2dec('8027b'));%���ִ���δ����ͨ�������ο����������ʹ�ã�����ֱ������sendָ��޻ض�ȷ�ϣ���
% dac.ReadValue(1,0);
% Attention:��P�ˣ�00000
classdef FTDAs <handle
    %��֪��<handle��ʲô���ã�ģ����û�С�
    properties
        ip; %�豸��ַ
        dac_handle;%�豸ָ��
        val_readout;
        str;%�����ַ���
        err_cnt;%ͨ�Ŵ����¼
%        nret;%���ݷ���ֵ
    end
    methods (Access = protected)
        function Open(obj)%ûë��
            fopen(obj.dac_handle);
        end
        
        function Close(obj)%ûë��
            fclose(obj.dac_handle);
        end   
    end
    
    methods
        function obj = FTDAs(ip)%Ҫ�޸����õĽ�����
            obj.dac_handle = tcpip(ip, 5000);
            set(obj.dac_handle,'Terminator','LF');%���ý�������Ӱ���ѯ�����������ò��ԣ����ܵ��¶�ȡ�����ʱ����Ŀǰ����ɶ������Ҳû���Գɹ���
            obj.err_cnt =0;
        end
        
        
        function result = ReadValue(obj,DA_id, DA_value)
            obj.Open();
            obj.str=sprintf('DA=%d;RW=0;ADDR=0x01;VAL=0x%05X',DA_id,DA_value);
            fprintf(obj.dac_handle, obj.str);
            obj.str= fscanf(obj.dac_handle,'%s',7);
            result =  sscanf(obj.str,'0x%05X'); %��obj.str���ַ���ת�������֣�Ŀǰ����ֵ�Ǹ�24λ16����������Ҫ�����ǰ4bit��֮����޸Ĺ̼���ֱ�ӷ���20λʮ������������0x����
            fprintf('Readback 0x%05X\n', result);%�����ã���ӡһ�»ض�ֵ
            obj.Close();
        end
        
        function send(obj, DA_id, DA_rw, DA_addr, DA_value)%ԭʼֱд����
            obj.Open();
            obj.str=sprintf('DA=%1d;RW=%1d;ADDR=0x%02X;VAL=0x%05X',DA_id,DA_rw, DA_addr, DA_value);
            fprintf(obj.dac_handle, obj.str);
            obj.Close();
        end
        
        function result = writeval(obj, DA_id, DA_value)%ֻд0x01��ַ��ֵ�������һ�λض������ض�ֵ�����ϲ����ý��жԱȼ��飬��С���ʻ�дʧ�ܡ�
            obj.Open();
            obj.str=sprintf('DA=%d;RW=1;ADDR=0x01;VAL=0x%05X',DA_id,DA_value);
            fprintf(obj.dac_handle, obj.str);%дһ��
            obj.str=sprintf('DA=%d;RW=0;ADDR=0x01;VAL=0x%05X',DA_id,DA_value);
            fprintf(obj.dac_handle, obj.str);%��һ��
            obj.val_readout= fscanf(obj.dac_handle,'%s',7);
            obj.str = sprintf('0x%05X',DA_value);
            if (~isempty(obj.val_readout) && strcmp(obj.str, obj.val_readout))
                result=0;
            else
                result=1;
                obj.err_cnt = obj.err_cnt+1;
            end
            obj.Close();
        end

        function SetValue(obj,DA_id, Value)
            while(writeval(obj,DA_id, Value))%whileѭ��ȷ��д�ɹ�������Ҳ�п�����Ϊ����ɳ���ѭ�������Կ�����ೢ��3�Σ������ʧ���򱨴��ء�
                fprintf('SetValule error occured, retrying...\n');
            end
            pause(0.2);
        end        
        
    end
end

