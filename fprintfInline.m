function formatSpec = fprintfInline(txtLabel,nIterations)
  
% build formatSpec: put necessary backspaces in, define format specifier (leftalign)

nBlanks = floor(log10(nIterations))+1;
fprintf([txtLabel blanks(nBlanks)]);
backString = '';

for n = 1:nBlanks
  backString = strcat(backString,'\b');
end
formatSpec = [backString '%-' num2str(nBlanks) 'd'];

end
