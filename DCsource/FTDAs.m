% author:F.Liang
% data:2018/4/1
% version:1.0
% filename:FTDAs.m
% describe:梁福田定制款高精度DC源操作类，封装版，每次读写进行开关tcp连接。
% 目前DA工作在tcpip server模式，接受的指令形式只有一种，字符串类型'DA=1;RW=1;ADDR=0x01;VAL=0x00000;'
% DA=(1,2,3,4);RW=(1W,0R);ADDR=0x(01,02,03,04);VAL=0x(00000~FFFFF);
% DA有4个通道，设备上每个通道对应一对互补反向输出，原则上正端输出的稳定性会较负端输出好，因为负端有外部网络电阻独立元件，可能影响性能。
% DA在接受RW=1;的写操作时，没有任何返回值，在接受RW=0;时，直接返回从DAC指定通道都回来的电压设置值。返回数据形式为字符串'0x00000\n'?结束符可后期调整硬件修改以配合软件设计。
% RW=0时，VAL无意义，但还是按规范格式带进去，可以传递0进入。
% ADDR常规用户只使用01，其他地址为DAC其他功能设置型寄存器，只对高级用户使用，所以常规读写函数不需要用户提供地址信息，函数内写死。
% VAL为电压设置值，20bit二进制码，0x00000最小，对应输出范围有可能是-2V~+2V，-3.5V~+3.5V，-7V~+7V，视电路直接配置情况（不能随意调节，需要开密封盖），目前版本为-2V~+2V
% 部分设备可能采用18bitDAC，设置值形式不变，传20bit到DAC，只是最低两bit无意义。目前版本为20bit。
% 使用范例：
% dac_ip='10.0.200.1';
% dac = FTDAs(dac_ip);
% dac.SetValue(1,hex2dec('8027b'));%部分代码未调试通，仅供参考，如果急需使用，建议直接运行send指令（无回读确认）。
% dac.ReadValue(1,0);
% Attention:接P端：00000
classdef FTDAs <handle
    %不知道<handle是什么作用，模板里没有。
    properties
        ip; %设备地址
        dac_handle;%设备指针
        val_readout;
        str;%传递字符串
        err_cnt;%通信错误记录
%        nret;%传递返回值
    end
    methods (Access = protected)
        function Open(obj)%没毛病
            fopen(obj.dac_handle);
        end
        
        function Close(obj)%没毛病
            fclose(obj.dac_handle);
        end   
    end
    
    methods
        function obj = FTDAs(ip)%要修改设置的结束符
            obj.dac_handle = tcpip(ip, 5000);
            set(obj.dac_handle,'Terminator','LF');%设置结束符，影响查询等命令，如果设置不对，可能导致读取结果超时，但目前设置啥结束符也没测试成功。
            obj.err_cnt =0;
        end
        
        
        function result = ReadValue(obj,DA_id, DA_value)
            obj.Open();
            obj.str=sprintf('DA=%d;RW=0;ADDR=0x01;VAL=0x%05X',DA_id,DA_value);
            fprintf(obj.dac_handle, obj.str);
            obj.str= fscanf(obj.dac_handle,'%s',7);
            result =  sscanf(obj.str,'0x%05X'); %将obj.str的字符串转换成数字，目前返回值是个24位16进制数，需要处理掉前4bit，之后会修改固件，直接返回20位十六进制数（带0x）。
            fprintf('Readback 0x%05X\n', result);%调试用，打印一下回读值
            obj.Close();
        end
        
        function send(obj, DA_id, DA_rw, DA_addr, DA_value)%原始直写函数
            obj.Open();
            obj.str=sprintf('DA=%1d;RW=%1d;ADDR=0x%02X;VAL=0x%05X',DA_id,DA_rw, DA_addr, DA_value);
            fprintf(obj.dac_handle, obj.str);
            obj.Close();
        end
        
        function result = writeval(obj, DA_id, DA_value)%只写0x01地址的值，并完成一次回读，将回读值用于上层设置进行对比检验，很小概率会写失败。
            obj.Open();
            obj.str=sprintf('DA=%d;RW=1;ADDR=0x01;VAL=0x%05X',DA_id,DA_value);
            fprintf(obj.dac_handle, obj.str);%写一次
            obj.str=sprintf('DA=%d;RW=0;ADDR=0x01;VAL=0x%05X',DA_id,DA_value);
            fprintf(obj.dac_handle, obj.str);%读一次
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
            while(writeval(obj,DA_id, Value))%while循环确保写成功，但是也有可能因为意外干成死循环，可以考虑最多尝试3次，如果都失败则报错返回。
                fprintf('SetValule error occured, retrying...\n');
            end
            pause(0.2);
        end        
        
    end
end

