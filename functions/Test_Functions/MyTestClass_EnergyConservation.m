%Testclass to check for energy conservation
classdef MyTestClass_EnergyConservation < matlab.unittest.TestCase
   properties
        t,f,fund_If,fund_It,shg_If,shg_It,spm_If,spm_It,xpm_If,xpm_It,xpm_If2,xpm_It2
        rel_tolerance=0.01;
        E_scal;
   end

    methods (Test)
        %%Tests for Energy Conservation   
        %% Tests in Frequency
        %Checks energy conservation between input pulse and propagated output pulse
        function test_EnergyConservationInOut_SPM(testCase)
            actSolution = get_Fluence(testCase.f,testCase.fund_If);
            expSolution = get_Fluence(testCase.f,testCase.spm_If);
            testCase.verifyEqual(round(actSolution,1),testCase.E_scal*round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end        
        %Checks energy conservation between input pulse and propagated output pulse
        function test_EnergyConservationInOutXPM_fund(testCase)
            actSolution = get_Fluence(testCase.f,testCase.fund_If);
            expSolution = get_Fluence(testCase.f,testCase.xpm_If);
            testCase.verifyEqual(round(actSolution,1),testCase.E_scal*round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end
        %Checks energy conservation between input pulse and propagated output pulse
        function test_EnergyConservationInOutXPM_SHG(testCase)
            actSolution = get_Fluence(testCase.f,testCase.shg_If);
            expSolution = get_Fluence(testCase.f,testCase.xpm_If2);
            testCase.verifyEqual(round(actSolution,1),testCase.E_scal*round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end
        
        %% Tests in Time
        %Checks energy conservation between input pulse and propagated output pulse
        function test_EnergyConservationInOut_SPM_t(testCase)
            actSolution = get_Fluence(testCase.t,testCase.fund_It);
            expSolution = get_Fluence(testCase.t,testCase.spm_It);
            testCase.verifyEqual(round(actSolution,1),testCase.E_scal*round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end        
        %Checks energy conservation between input pulse and propagated output pulse
        function test_EnergyConservationInOutXPM_fund_t(testCase)
            actSolution = get_Fluence(testCase.t,testCase.fund_It);
            expSolution = get_Fluence(testCase.t,testCase.xpm_It);
            testCase.verifyEqual(round(actSolution,1),testCase.E_scal*round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end
        %Checks energy conservation between input pulse and propagated output pulse
        function test_EnergyConservationInOutXPM_SHG_t(testCase)
            actSolution = get_Fluence(testCase.t,testCase.shg_It);
            expSolution = get_Fluence(testCase.t,testCase.xpm_It2);
            testCase.verifyEqual(round(actSolution,1),testCase.E_scal*round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end        
                
        %% Tests with respect to Fourier Transform t<->f  
        %% In
        %Checks energy conservation of fund pulse 
        function test_EnergyConservationIn_fund(testCase)
            actSolution = get_Fluence(testCase.f,testCase.fund_If);
            expSolution = get_Fluence(testCase.t,testCase.fund_It);
            testCase.verifyEqual(round(actSolution,1),round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end 
        %Checks energy conservation of shg pulse 
        function test_EnergyConservationIn_SHG(testCase)
            actSolution = get_Fluence(testCase.f,testCase.shg_If);
            expSolution = get_Fluence(testCase.t,testCase.shg_It);
            testCase.verifyEqual(round(actSolution,1),round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end  
        
        %% Out
        %Checks energy conservation of spm pulse
        function test_EnergyConservationOut_SPM(testCase)
            actSolution = get_Fluence(testCase.f,testCase.spm_If);
            expSolution = get_Fluence(testCase.t,testCase.spm_It);
            testCase.verifyEqual(round(actSolution,1),round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end
        %Checks energy conservation of xpm fundamental pulse
        function test_EnergyConservationOut_XPM_fund(testCase)
            actSolution = get_Fluence(testCase.f,testCase.xpm_If);
            expSolution = get_Fluence(testCase.t,testCase.xpm_It);
            testCase.verifyEqual(round(actSolution,1),round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end
        %Checks energy conservation of xpm SHG pulse
        function test_EnergyConservationOut_XPM_SHG(testCase)
            actSolution = get_Fluence(testCase.f,testCase.xpm_If2);
            expSolution = get_Fluence(testCase.t,testCase.xpm_It2);
            testCase.verifyEqual(round(actSolution,1),round(expSolution,1),'RelTol',testCase.rel_tolerance);
        end        
        
    end
    
end 


