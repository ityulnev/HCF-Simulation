%Constructor for testCase Energy cosnervation in HCF simulation
function answer=test_energy(varargin)

testCase = MyTestClass_EnergyConservation;

testCase.t=varargin{1}; 
testCase.f=varargin{2}; 
testCase.fund_If=varargin{3}; 
testCase.fund_It=varargin{4}; 
testCase.shg_If=varargin{5};
testCase.shg_It=varargin{6}; 
testCase.spm_If=varargin{7}; 
testCase.spm_It=varargin{8}; 
testCase.xpm_If=varargin{9}; 
testCase.xpm_It=varargin{10}; 
testCase.xpm_If2=varargin{11}; 
testCase.xpm_It2=varargin{12}; 
testCase.E_scal=varargin{13}; 

answer = run(testCase);
end