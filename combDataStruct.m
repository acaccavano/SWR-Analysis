function data = combDataStruct(data1, data2, data3)
if (nargin < 3) data3 = []; end
if (nargin < 2) data2 = []; end
if (nargin < 1) data1 = []; end

% Initialize data structure and field names
if ~isempty(data1)
  f1 = fieldnames(data1);
  data(1) = data1;
  
  % combining data2 structure:
  if ~isempty(data2)
    f2 = fieldnames(data2);
    diff = setdiff(f1, f2);
    if ~isempty(diff)
      for i = 1:length(diff)
        if ~isfield(data(1), diff{i}) data(1).(diff{i}) = struct; end
        if ~isfield(data2, diff{i}) data2.(diff{i}) = struct; end
      end
    end
    data(2) = data2;
    
    % combining data3 structure:
    if ~isempty(data3)
      f1 = fieldnames(data(1));
      f3 = fieldnames(data3);
      diff = setdiff(f1, f3);
      if ~isempty(diff)
        for i = 1:length(diff)
          if ~isfield(data(1), diff{i}) data(1).(diff{i}) = struct; end
          if ~isfield(data(2), diff{i}) data(2).(diff{i}) = struct; end
          if ~isfield(data3, diff{i}) data3.(diff{i}) = struct; end
        end
      end
      data(3) = data3;
      
    end
  end
end

end
