boots = 100
threshold = 0.1
parpool(9)


addpath(genpath('.'))

fileNamesList = {'filtered_ibd_t_t_alignment_sr7d','filtered_ibd_t_alignment_sr7d','filtered_ibd_t_g_alignment_sr7d','filtered_ibd_t_m_alignment_sr7d','filtered_ibd_t_g_m_alignment_sr7d','filtered_ibd_t_t_noalignment_sr7d','filtered_ibd_t_noalignment_sr7d','filtered_ibd_t_g_noalignment_sr7d','filtered_ibd_t_m_noalignment_sr7d','filtered_ibd_t_g_m_noalignment_sr7d'}


bnpathscript 
for maxParents = 3:3
	parfor i = 1:8
		fileName = fileNamesList{i}
	    if contains(fileName,'_t_g_m')
	        allowedMatrix = 'HEGTM_Skeleton'
	    elseif  contains(fileName,'_t_g')
	        allowedMatrix = 'TG'
	    elseif  contains(fileName,'_t_m')
	        allowedMatrix = 'TM'
	    else
	        allowedMatrix = 'HEGTM'
	    end
		bootstrapLearnHeterogeneousDynamicBayesNetwork(fileName, '', boots, threshold,[], maxParents, false, false, true,1,1,1,0,4,false,allowedMatrix);  
	end
end



