% Read dataset
function [Tpower, Tx] = readData(file_tx, file_rss)
    % set file directory
%     file_tx : Tx坐标 
%     file_rss ：RSS信息

    % 导入选项
    Tx = readtable(file_tx, 'ReadVariableNames', false);
    Tx_num = height(Tx);
    opts = detectImportOptions(sprintf(file_rss,1), ...
        'FileType','text','CommentStyle','#','Delimiter',{' ','\t'});
    opts.VariableNames = {'Idx','X','Y','Z','Distance','Power','Phase'};
    opts.VariableTypes = repmat("double",1,7);

    % 读取首文件
    T0 = readtable(sprintf(file_rss,1), opts);
    P = zeros(height(T0), Tx_num);
    P(:,1) = T0.Power;

    % 读取剩余
    for k = 2:Tx_num
        Tt = readtable(sprintf(file_rss,k), opts);
        P(:,k) = Tt.Power;
    end

    % 构造输出 table
    T0.X = ceil(T0.X); 
    T0.Y = ceil(T0.Y); 
    T0.Z = ceil(T0.Z);
    Tpower = table(T0.Idx, T0.X, T0.Y, T0.Z, P, ...
        'VariableNames',{'Idx','X','Y','Z','Power'});
end
