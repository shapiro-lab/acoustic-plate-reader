clc

disp('Searching for Grid Files');
files = subdir('*Scan_Well*.mat');
disp(' ');

%%
clc


for i = 1:numel(files);
   fn = files(i).name;
   iname = [fn(1:end-3) 'png'];
   
   if exist(iname, 'file') ~= 2
   if contains(fn, 'Grid_SG')
       disp(fn)
       disp('On exclusion list');
       disp(' ');
   else
       
       
   disp(fn);
   disp('- Loading File...');
   load(fn);
   

   disp('- Analyzing...');
   
   try
       
    ax1 = (-(params.Scan.dim1_total - 1)/2:1:(params.Scan.dim1_total-1)/2) ...
    .* params.Scan.dim1_step .* params.Stages.step_distance * 1000;

    ax2 = (-(params.Scan.dim2_total - 1)/2:1:(params.Scan.dim2_total-1)/2) ...
    .* params.Scan.dim2_step .* params.Stages.step_distance * 1000;

    figure(1); clf;
    imagesc(ax1, ax2, reshape(params.Scan.Objective, params.Scan.dim2_total, params.Scan.dim1_total));
    xlabel('Distance (mm)');
    ylabel('Distance (mm)');
    axis image;
    colorbar;
    title([params.Name ' ' params.Time], 'Interpreter', 'None');
    
    iname = [fn(1:end-3) 'png'];
    print(1, '-dpng', iname);
    
    disp(' - Success!');    
   catch
       
       try
       
        disp(sprintf(' - Location Points: %1.0f', numel(params.Scan.Location(1,:))));
           
       % Find static dimension
       staticdim = 0;
       for j = 1:3
           if std(params.Scan.Location(j,:)) == 0
               staticdim = j;
           end
       end
       if staticdim == 0;
           throw(MException('','Could not find static dim'));
       end
       
       dim1 = 0;
       for j = 1:3
           if j ~= staticdim
               if params.Scan.Location(j,2) ~= params.Scan.Location(j,1)
                   dim1 = j;
               end
           end
       end
       
       if dim1 == 0;
           throw(MException('','Could not find dim1'));
       end
       
       for j = 1:3
           if j ~= dim1
               if j ~= staticdim
                   dim2 = j;
               end
           end
       end
       
       
       n1 = 0;
       for j = 1:numel(params.Scan.Location(1,:))-1;
           if n1 == 0
           if params.Scan.Location(dim2,j) ~= params.Scan.Location(dim2, j+1);
               n1 = j;
           end
           end
       end
       
       if n1 == 0;
           throw(MException('','Could not find n1'));
       end
       
       n2 = numel(params.Scan.Location(1,:)) / n1;
       
       if isfield(params.Scan, 'Objective');
           obj = params.Scan.Objective;
       elseif isfield(params.Scan, 'PeakToPeak');
           obj = params.Scan.PeakToPeak;
       elseif isfield(params.Scan, 'FFT_peak');
           obj = params.Scan.FFT_peak;
       elseif isfield(params.Result, 'Objective');
           obj = params.Results.Objective;
       else
           throw(MException('','Cannot find objective function'))
       end
       
       
              
       ds1 = params.Scan.Location(dim1,2     ) - params.Scan.Location(dim1,1);
       ds2 = params.Scan.Location(dim2,n1 + 1) - params.Scan.Location(dim2,1);
       
        ax1 = (-(n1 - 1)/2:1:(n1-1)/2) ...
        .* ds1 .* params.Stages.step_distance * 1000;

        ax2 = (-(n2 - 1)/2:1:(n2-1)/2) ...
        .* ds2 .* params.Stages.step_distance * 1000;

        figure(1); clf;
        imagesc(ax1, ax2, reshape(obj, n2, n1));
        xlabel('Distance (mm)');
        ylabel('Distance (mm)');
        axis image;
        colorbar;
        title([params.Name ' ' params.Time], 'Interpreter', 'None');

        figure(2); clf;
        scatter3(params.Scan.Location(1,:), params.Scan.Location(2,:), params.Scan.Location(3,:), 25*ones(numel(obj),1), obj);
        
        iname = [fn(1:end-3) 'png'];
        print(1, '-dpng', iname);

        pause(1);
        disp(' - Success!');  
        
       catch ex
           disp(' - Error')
           disp([' - ' ex.message]);
       end
        
   end
   
   disp(' ');
   
   end
   end
end