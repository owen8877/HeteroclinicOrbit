function sol = mergeInit()
    AB = load('backup/matlabAB.mat');
    BC = load('backup/matlabBC.mat');
    CD = load('backup/matlabCD.mat');
    sol = [AB.sol BC.sol CD.sol];
    sol = sol(:, 1:2:end);
end