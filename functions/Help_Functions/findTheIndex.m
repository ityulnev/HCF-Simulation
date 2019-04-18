%Find me the Index corresponding to the value
function index=findTheIndex(array,value)
[~,index] = min(abs( array-value ));
end