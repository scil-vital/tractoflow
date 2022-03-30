#!/usr/bin/env nextflow

import groovy.json.*

params.input = false
params.bids = false
params.bids_config = false
params.help = false
params.dti_shells = false
params.fodf_shells = false

if(params.help) {
    usage = file("$baseDir/USAGE")

    cpu_count = Runtime.runtime.availableProcessors()
    bindings = ["clean_bids":"$params.clean_bids",
                "b0_thr_extract_b0":"$params.b0_thr_extract_b0",
                "dwi_shell_tolerance":"$params.dwi_shell_tolerance",
                "run_resample_dwi":"$params.run_resample_dwi",
                "dwi_resolution":"$params.dwi_resolution",
                "dwi_interpolation":"$params.dwi_interpolation",
                "fa":"$params.fa",
                "min_fa":"$params.min_fa",
                "roi_radius":"$params.roi_radius",
                "set_frf":"$params.set_frf",
                "manual_frf":"$params.manual_frf",
                "mean_frf":"$params.mean_frf",
                "sh_order":"$params.sh_order",
                "basis":"$params.basis",
                "fodf_metrics_a_factor":"$params.fodf_metrics_a_factor",
                "relative_threshold":"$params.relative_threshold",
                "max_fa_in_ventricle":"$params.max_fa_in_ventricle",
                "min_md_in_ventricle":"$params.min_md_in_ventricle",
                "run_pft_tracking":"$params.run_pft_tracking",
                "pft_seeding_mask_type":"$params.pft_seeding_mask_type",
                "pft_fa_seeding_mask_threshold":"$params.pft_fa_seeding_mask_threshold",
                "pft_algo":"$params.pft_algo",
                "pft_seeding":"$params.pft_seeding",
                "pft_nbr_seeds":"$params.pft_nbr_seeds",
                "pft_step":"$params.pft_step",
                "pft_theta":"$params.pft_theta",
                "pft_min_len":"$params.pft_min_len",
                "pft_max_len":"$params.pft_max_len",
                "pft_compress_streamlines":"$params.pft_compress_streamlines",
                "pft_compress_value":"$params.pft_compress_value",
                "local_seeding_mask_type":"$params.local_seeding_mask_type",
                "local_fa_seeding_mask_threshold":"$params.local_fa_seeding_mask_threshold",
                "local_tracking_mask_type":"$params.local_tracking_mask_type",
                "local_fa_tracking_mask_threshold":"$params.local_fa_tracking_mask_threshold",
                "run_local_tracking":"$params.run_local_tracking",
                "local_compress_streamlines":"$params.local_compress_streamlines",
                "pft_random_seed":"$params.pft_random_seed",
                "local_algo":"$params.local_algo",
                "local_seeding":"$params.local_seeding",
                "local_nbr_seeds":"$params.local_nbr_seeds",
                "local_step":"$params.local_step",
                "local_theta":"$params.local_theta",
                "local_sfthres":"$params.local_sfthres",
                "local_sfthres_init":"$params.local_sfthres_init",
                "local_min_len":"$params.local_min_len",
                "local_max_len":"$params.local_max_len",
                "local_compress_value":"$params.local_compress_value",
                "local_random_seed":"$params.local_random_seed",
                "cpu_count":"$cpu_count",
                "processes_fodf":"$params.processes_fodf"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "TractoFlow-minimal pipeline"
log.info "==================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}


if (params.input && !(params.bids && params.bids_config)){
    log.info "Input: $params.input"
    root = file(params.input)
    data = Channel
        .fromFilePairs("$root/**/*{bval,bvec,dwi.nii.gz,mask_wm.nii.gz}",
                       size: 4,
                       maxDepth:1,
                       flat: true) {it.parent.name}

    data
        .into{in_data; check_subjects_number}

    pft_maps = Channel
        .fromFilePairs("$root/**/*{map_csf.nii.gz,map_gm.nii.gz,map_wm.nii.gz}",
                       size: 3,
                       maxDepth:1,
                       flat: true) {it.parent.name}

}
else if (params.bids || params.bids_config){
    if (!params.bids_config) {
        log.info "Input BIDS: $params.bids"
        if (params.participants_label) {
            Integer start = workflow.commandLine.indexOf("participants_label") + "participants_label".length();
            participant_cleaned = workflow.commandLine.substring(start, workflow.commandLine.indexOf("--", start) == -1 ? workflow.commandLine.length() : workflow.commandLine.indexOf("--", start)).replace("=", "").replace("\'", "")
            log.info "Participants: $participant_cleaned"
        }
        log.info "Clean_bids: $params.clean_bids"
        log.info ""

        bids = file(params.bids)

        process Read_BIDS {
            publishDir = params.Read_BIDS_Publish_Dir
            scratch = false
            stageInMode = 'symlink'
            tag = {"Read_BIDS"}
            errorStrategy = { task.attempt <= 3 ? 'retry' : 'terminate' }

            input:
            file(bids_folder) from bids

            output:
            file "tractoflow_bids_struct.json" into bids_struct

            script:
            participants_flag =\
            params.participants_label ? '--participants_label ' + participant_cleaned : ""

            clean_flag = params.clean_bids ? '--clean ' : ''

            """
            scil_validate_bids.py $bids_folder tractoflow_bids_struct.json\
                --readout $params.readout $participants_flag $clean_flag
            """
        }
    }

    else {
        log.info "BIDS config: $params.bids_config"
        config = file(params.bids_config)
        bids_struct = Channel.from(config)
    }

    ch_in_data = Channel.create()
    bids_struct.map{it ->
    jsonSlurper = new JsonSlurper()
        data = jsonSlurper.parseText(it.getText())
        for (item in data){
            sid = "sub-" + item.subject

            if (item.session){
                sid += "_ses-" + item.session
            }

            if (item.run){
                sid += "_run-" + item.run
            }
            for (key in item.keySet()){
                if(item[key] == 'todo'){
                    error "Error ~ Please look at your tractoflow_bids_struct.json " +
                    "in Read_BIDS folder.\nPlease fix todo fields and give " +
                    "this file in input using --bids_config option instead of" +
                    "using --bids."
                }
                else if (item[key] == 'error_readout'){
                    error "Error ~ Please look at your tractoflow_bids_struct.json " +
                    "in Read_BIDS folder.\nPlease fix error_readout fields. "+
                    "This error indicate that readout time looks wrong.\n"+
                    "Please correct the value or remove the subject in the json and " +
                    "give the updated file in input using --bids_config option instead of" +
                    "using --bids."
                }
            }
            sub = [sid, file(item.bval), file(item.bvec), file(item.dwi),
                   file(item.t1), item.TotalReadoutTime, item.DWIPhaseEncodingDir[0]]
            ch_in_data.bind(sub)
        }

        ch_in_data.close()
    }

    ch_in_data.into{in_data; check_subjects_number}
}
else {
    error "Error ~ Please use --input, --bids or --bids_config for the input data."
}

if (!params.dti_shells || !params.fodf_shells){
    error "Error ~ Please set the DTI and fODF shells to use."
}


if (params.pft_seeding_mask_type != "wm" && params.pft_seeding_mask_type != "interface" && params.pft_seeding_mask_type != "fa"){
    error "Error ~ --pft_seeding_mask_type can only take wm, interface or fa. Please select one of these choices"
}

if (params.local_seeding_mask_type != "wm" && params.local_seeding_mask_type != "fa"){
    error "Error ~ --local_seeding_mask_type can only take wm or fa. Please select one of these choices"
}

if (params.local_tracking_mask_type != "wm" && params.local_tracking_mask_type != "fa"){
    error "Error ~ --local_tracking_mask_type can only take wm or fa. Please select one of these choices"
}

if (params.local_algo != "det" && params.local_algo != "prob"){
    error "Error ~ --local_algo can only take det or prob. Please select one of these choices"
}

if (params.pft_algo != "det" && params.pft_algo != "prob"){
    error "Error ~ --pft_algo can only take det or prob. Please select one of these choices"
}

if (params.local_seeding != "nt" && params.local_seeding != "npv"){
    error "Error ~ --local_seeding can only take nt or npv. Please select one of these choices"
}

if (params.pft_seeding != "nt" && params.pft_seeding != "npv"){
    error "Error ~ --pft_seeding can only take nt or npv. Please select one of these choices"
}

if (params.run_pft_tracking && workflow.profile.contains("ABS")){
    error "Error ~ PFT tracking cannot be run with Atlas Based Segmentation (ABS) profile"
}


(dwi_mask_for_resample, gradients, dwi_gradients_for_extract_b0) = in_data
    .map{sid, bvals, bvecs, dwi, mask_wm -> [
                                        tuple(sid, dwi, mask_wm),
                                        tuple(sid, bvals, bvecs),
                                        tuple(sid, dwi, bvals, bvecs)]}
    .separate(3)

gradients.into{gradients_for_dti_shell; gradients_for_fodf_shell; gradients_for_extract_b0}

check_subjects_number.count().set{ number_subj_for_null_check }

number_subj_for_null_check
.subscribe{a -> if (a == 0)
    error "Error ~ No subjects found. Please check the naming convention, your --input path or your BIDS folder."}

if (params.set_frf && params.mean_frf){
    error "Error ~ --set_frf and --mean_frf are activated. Please choose only one of these options. "
}

if (params.pft_random_seed instanceof String){
    pft_random_seed = params.pft_random_seed?.tokenize(',')
}
else{
    pft_random_seed = params.pft_random_seed
}

if (params.local_random_seed instanceof String){
    local_random_seed = params.local_random_seed?.tokenize(',')
}
else{
    local_random_seed = params.local_random_seed
}


process README {
    cpus 1
    publishDir = params.Readme_Publish_Dir
    tag = "README"

    output:
    file "readme.txt"

    script:
    String list_options = new String();
    for (String item : params) {
        list_options += item + "\n"
    }
    """
    echo "TractoFlow pipeline\n" >> readme.txt
    echo "Start time: $workflow.start\n" >> readme.txt
    echo "[Command-line]\n$workflow.commandLine\n" >> readme.txt
    echo "[Git Info]\n" >> readme.txt
    echo "$workflow.repository - $workflow.revision [$workflow.commitId]\n" >> readme.txt
    echo "[Options]\n" >> readme.txt
    echo "$list_options" >> readme.txt
    """
}

process Resample_DWI_And_WM_Mask {
    cpus 3

    input:
    set sid, file(dwi), file(wm_mask) from dwi_mask_for_resample

    output:
    set sid, "${sid}__dwi_resampled.nii.gz" into\
        dwi_for_extract_b0,
        dwi_for_extract_dti_shell,
        dwi_for_extract_fodf_shell
    set sid, "${sid}__mask_wm_resampled.nii.gz" into wm_mask

    script:
    if (params.run_resample_dwi)
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_resample_volume.py $dwi \
            dwi_resample.nii.gz \
            --resolution $params.dwi_resolution \
            --interp  $params.dwi_interpolation
        fslmaths dwi_resample.nii.gz -thr 0 dwi_resample_clipped.nii.gz
        mv dwi_resample_clipped.nii.gz ${sid}__dwi_resampled.nii.gz

        scil_resample_volume.py $wm_mask "${sid}__mask_wm_resampled.nii.gz" \
        --resolution $params.dwi_resolution --interp nn
        """
    else
        """
        mv $dwi ${sid}__dwi_resampled.nii.gz
        mv $wm_mask "${sid}__mask_wm_resampled.nii.gz"
        """
}

process Resample_PFT_Maps {
    cpus 3

    input:
    set sid, file(csf), file(gm), file(wm) from pft_maps

    output:
    set sid,
    "${sid}__map_csf_resampled.nii.gz",
    "${sid}__map_gm_resampled.nii.gz",
     "${sid}__map_wm_resampled.nii.gz" into pft_maps_resampled

    script:
    if (params.run_resample_dwi)
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1

        scil_resample_volume.py $csf "${sid}__map_csf_resampled.nii.gz" \
        --resolution $params.dwi_resolution --interp lin
        scil_resample_volume.py $gm "${sid}__map_gm_resampled.nii.gz" \
        --resolution $params.dwi_resolution --interp lin
        scil_resample_volume.py $wm "${sid}__map_wm_resampled.nii.gz" \
        --resolution $params.dwi_resolution --interp lin
        """
    else
        """
        mv $csf ${sid}__map_csf_resampled.nii.gz
        mv $gm ${sid}__map_gm_resampled.nii.gz
        mv $wm ${sid}__map_wm_resampled.nii.gz
        """
}

dwi_for_extract_b0
    .join(gradients_for_extract_b0)
    .set{dwi_and_grad_for_extract_b0}

process Compute_B0_Mask {
    cpus 3

    input:
    set sid, file(dwi), file(bval), file(bvec) from dwi_and_grad_for_extract_b0

    output:
    set sid, "${sid}__b0_mask_resampled.nii.gz" into\
        b0_mask_for_dti_metrics,
        b0_mask_for_fodf,
        b0_mask_for_rf

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0_resampled.nii.gz --mean\
        --b0_thr $params.b0_thr_extract_b0 --force_b0_threshold
    mrthreshold ${sid}__b0_resampled.nii.gz ${sid}__b0_mask_resampled.nii.gz\
        --abs 0.00001 -nthreads 1
    """
}


dwi_for_extract_dti_shell
    .join(gradients_for_dti_shell)
    .set{dwi_and_grad_for_extract_dti_shell}

process Extract_DTI_Shell {
    cpus 3

    input:
    set sid, file(dwi), file(bval), file(bvec)\
        from dwi_and_grad_for_extract_dti_shell

    output:
    set sid, "${sid}__dwi_dti.nii.gz", "${sid}__bval_dti",
        "${sid}__bvec_dti" into \
        dwi_and_grad_for_dti_metrics, \
        dwi_and_grad_for_rf

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_extract_dwi_shell.py $dwi \
        $bval $bvec $params.dti_shells ${sid}__dwi_dti.nii.gz \
        ${sid}__bval_dti ${sid}__bvec_dti -t $params.dwi_shell_tolerance -f
    """
}

dwi_and_grad_for_dti_metrics
    .join(b0_mask_for_dti_metrics)
    .set{dwi_and_grad_for_dti_metrics}

process DTI_Metrics {
    cpus 3

    input:
    set sid, file(dwi), file(bval), file(bvec), file(b0_mask)\
        from dwi_and_grad_for_dti_metrics

    output:
    file "${sid}__ad.nii.gz"
    file "${sid}__evecs.nii.gz"
    file "${sid}__evecs_v1.nii.gz"
    file "${sid}__evecs_v2.nii.gz"
    file "${sid}__evecs_v3.nii.gz"
    file "${sid}__evals.nii.gz"
    file "${sid}__evals_e1.nii.gz"
    file "${sid}__evals_e2.nii.gz"
    file "${sid}__evals_e3.nii.gz"
    file "${sid}__fa.nii.gz"
    file "${sid}__ga.nii.gz"
    file "${sid}__rgb.nii.gz"
    file "${sid}__md.nii.gz"
    file "${sid}__mode.nii.gz"
    file "${sid}__norm.nii.gz"
    file "${sid}__rd.nii.gz"
    file "${sid}__tensor.nii.gz"
    file "${sid}__nonphysical.nii.gz"
    file "${sid}__pulsation_std_dwi.nii.gz"
    file "${sid}__residual.nii.gz"
    file "${sid}__residual_iqr_residuals.npy"
    file "${sid}__residual_mean_residuals.npy"
    file "${sid}__residual_q1_residuals.npy"
    file "${sid}__residual_q3_residuals.npy"
    file "${sid}__residual_residuals_stats.png"
    file "${sid}__residual_std_residuals.npy"
    set sid, "${sid}__fa.nii.gz", "${sid}__md.nii.gz" into fa_md_for_fodf
    set sid, "${sid}__fa.nii.gz" into\
        fa_for_pft_tracking, fa_for_local_tracking_mask, fa_for_local_seeding_mask

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_dti_metrics.py $dwi $bval $bvec --mask $b0_mask\
        --ad ${sid}__ad.nii.gz --evecs ${sid}__evecs.nii.gz\
        --evals ${sid}__evals.nii.gz --fa ${sid}__fa.nii.gz\
        --ga ${sid}__ga.nii.gz --rgb ${sid}__rgb.nii.gz\
        --md ${sid}__md.nii.gz --mode ${sid}__mode.nii.gz\
        --norm ${sid}__norm.nii.gz --rd ${sid}__rd.nii.gz\
        --tensor ${sid}__tensor.nii.gz\
        --non-physical ${sid}__nonphysical.nii.gz\
        --pulsation ${sid}__pulsation.nii.gz\
        --residual ${sid}__residual.nii.gz\
        -f --force_b0_threshold
    """
}

dwi_for_extract_fodf_shell
    .join(gradients_for_fodf_shell)
    .set{dwi_and_grad_for_extract_fodf_shell}

process Extract_FODF_Shell {
    cpus 3

    input:
    set sid, file(dwi), file(bval), file(bvec)\
        from dwi_and_grad_for_extract_fodf_shell

    output:
    set sid, "${sid}__dwi_fodf.nii.gz", "${sid}__bval_fodf",
        "${sid}__bvec_fodf" into\
        dwi_and_grad_for_fodf

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_extract_dwi_shell.py $dwi \
        $bval $bvec $params.fodf_shells ${sid}__dwi_fodf.nii.gz \
        ${sid}__bval_fodf ${sid}__bvec_fodf -t $params.dwi_shell_tolerance -f
    """
}

wm_mask
    .into{wm_mask_for_pft_tracking;wm_mask_for_local_tracking_mask;wm_mask_for_local_seeding_mask}

dwi_and_grad_for_rf
    .join(b0_mask_for_rf)
    .set{dwi_b0_for_rf}

process Compute_FRF {
    cpus 3

    input:
    set sid, file(dwi), file(bval), file(bvec), file(b0_mask)\
        from dwi_b0_for_rf

    output:
    set sid, "${sid}__frf.txt" into unique_frf, unique_frf_for_mean
    file "${sid}__frf.txt" into all_frf_to_collect

    script:
    if (params.set_frf)
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_ssst_frf.py $dwi $bval $bvec frf.txt --mask $b0_mask\
        --fa $params.fa --min_fa $params.min_fa --min_nvox $params.min_nvox\
        --roi_radii $params.roi_radius --force_b0_threshold
        scil_set_response_function.py frf.txt $params.manual_frf ${sid}__frf.txt
        """
    else
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_ssst_frf.py $dwi $bval $bvec ${sid}__frf.txt --mask $b0_mask\
        --fa $params.fa --min_fa $params.min_fa --min_nvox $params.min_nvox\
        --roi_radii $params.roi_radius --force_b0_threshold
        """
}

all_frf_to_collect
    .collect()
    .set{all_frf_for_mean_frf}

process Mean_FRF {
    cpus 1
    publishDir = params.Mean_FRF_Publish_Dir
    tag = {"All_FRF"}

    input:
    file(all_frf) from all_frf_for_mean_frf

    output:
    file "mean_frf.txt" into mean_frf

    when:
    params.mean_frf && !params.set_frf

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_mean_frf.py $all_frf mean_frf.txt
    """
}

frf_for_fodf = unique_frf

if (params.mean_frf) {
    frf_for_fodf = unique_frf_for_mean
                   .merge(mean_frf)
                   .map{it -> [it[0], it[2]]}
}

dwi_and_grad_for_fodf
    .join(b0_mask_for_fodf)
    .join(fa_md_for_fodf)
    .join(frf_for_fodf)
    .set{dwi_b0_metrics_frf_for_fodf}

process FODF_Metrics {
    cpus params.processes_fodf

    input:
    set sid, file(dwi), file(bval), file(bvec), file(b0_mask), file(fa),
        file(md), file(frf) from dwi_b0_metrics_frf_for_fodf

    output:
    set sid, "${sid}__fodf.nii.gz" into fodf_for_pft_tracking, fodf_for_local_tracking
    file "${sid}__peaks.nii.gz"
    file "${sid}__peak_indices.nii.gz"
    file "${sid}__afd_max.nii.gz"
    file "${sid}__afd_total.nii.gz"
    file "${sid}__afd_sum.nii.gz"
    file "${sid}__nufo.nii.gz"

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_ssst_fodf.py $dwi $bval $bvec $frf ${sid}__fodf.nii.gz\
        --sh_order $params.sh_order --sh_basis $params.basis --force_b0_threshold\
        --mask $b0_mask --processes $task.cpus

    scil_compute_fodf_max_in_ventricles.py ${sid}__fodf.nii.gz $fa $md\
        --max_value_output ventricles_fodf_max_value.txt --sh_basis $params.basis\
        --fa_t $params.max_fa_in_ventricle --md_t $params.min_md_in_ventricle\
        -f

    a_threshold=\$(echo $params.fodf_metrics_a_factor*\$(cat ventricles_fodf_max_value.txt)|bc)

    scil_compute_fodf_metrics.py ${sid}__fodf.nii.gz\
        --mask $b0_mask --sh_basis $params.basis\
        --peaks ${sid}__peaks.nii.gz --peak_indices ${sid}__peak_indices.nii.gz\
        --afd_max ${sid}__afd_max.nii.gz --afd_total ${sid}__afd_total.nii.gz\
        --afd_sum ${sid}__afd_sum.nii.gz --nufo ${sid}__nufo.nii.gz\
        --rt $params.relative_threshold --at \${a_threshold}
    """
}

pft_maps_resampled
    .map{ sid, csf, gm, wm -> [sid, wm, gm, csf]}
    .set{ map_wm_gm_csf_for_pft_maps }

process PFT_Tracking_Maps {
    cpus 1

    input:
    set sid, file(wm), file(gm), file(csf) from map_wm_gm_csf_for_pft_maps

    output:
    set sid, "${sid}__map_include.nii.gz",
        "${sid}__map_exclude.nii.gz" into pft_maps_for_pft_tracking
    set sid, "${sid}__interface.nii.gz" into interface_for_pft_seeding_mask

    when:
        params.run_pft_tracking

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_maps_for_particle_filter_tracking.py $wm $gm $csf \
        --include ${sid}__map_include.nii.gz \
        --exclude ${sid}__map_exclude.nii.gz \
        --interface ${sid}__interface.nii.gz -f
    """
}
wm_mask_for_pft_tracking
    .join(fa_for_pft_tracking)
    .join(interface_for_pft_seeding_mask)
    .set{wm_fa_int_for_pft}

process PFT_Seeding_Mask {
    cpus 1

    input:
    set sid, file(wm), file(fa), file(interface_mask) from wm_fa_int_for_pft

    output:
    set sid, "${sid}__pft_seeding_mask.nii.gz" into seeding_mask_for_pft

    when:
        params.run_pft_tracking

    script:
    if (params.pft_seeding_mask_type == "wm")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_image_math.py union $wm $interface_mask ${sid}__pft_seeding_mask.nii.gz\
            --data_type uint8
        """
    else if (params.pft_seeding_mask_type == "interface")
        """
        mv $interface_mask ${sid}__pft_seeding_mask.nii.gz
        """
    else if (params.pft_seeding_mask_type == "fa")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrcalc $fa $params.pft_fa_seeding_mask_threshold -ge ${sid}__pft_seeding_mask.nii.gz\
          -datatype uint8
        """
}

fodf_for_pft_tracking
    .join(pft_maps_for_pft_tracking)
    .join(seeding_mask_for_pft)
    .set{fodf_maps_for_pft_tracking}

process PFT_Tracking {
    cpus 2

    input:
    set sid, file(fodf), file(include), file(exclude), file(seed)\
        from fodf_maps_for_pft_tracking
    each curr_seed from pft_random_seed

    output:
    file "${sid}__pft_tracking_${params.pft_algo}_${params.pft_seeding_mask_type}_seed_${curr_seed}.trk"

    when:
        params.run_pft_tracking

    script:
    compress =\
        params.pft_compress_streamlines ? '--compress ' + params.pft_compress_value : ''
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_pft.py $fodf $seed $include $exclude\
            ${sid}__pft_tracking_${params.pft_algo}_${params.pft_seeding_mask_type}_seed_${curr_seed}.trk\
            --algo $params.pft_algo --$params.pft_seeding $params.pft_nbr_seeds\
            --seed $curr_seed --step $params.pft_step --theta $params.pft_theta\
            --sfthres $params.pft_sfthres --sfthres_init $params.pft_sfthres_init\
            --min_length $params.pft_min_len --max_length $params.pft_max_len\
            --particles $params.pft_particles --back $params.pft_back\
            --forward $params.pft_front $compress --sh_basis $params.basis
        """
}

wm_mask_for_local_tracking_mask
    .join(fa_for_local_tracking_mask)
    .set{wm_fa_for_local_tracking_mask}

process Local_Tracking_Mask {
    cpus 1

    input:
    set sid, file(wm), file(fa) from wm_fa_for_local_tracking_mask

    output:
    set sid, "${sid}__local_tracking_mask.nii.gz" into tracking_mask_for_local

    when:
        params.run_local_tracking

    script:
    if (params.local_tracking_mask_type == "wm")
        """
        mv $wm ${sid}__local_tracking_mask.nii.gz
        """
    else if (params.local_tracking_mask_type == "fa")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrcalc $fa $params.local_fa_tracking_mask_threshold -ge ${sid}__local_tracking_mask.nii.gz\
          -datatype uint8
        """
}

wm_mask_for_local_seeding_mask
    .join(fa_for_local_seeding_mask)
    .set{wm_fa_for_local_seeding_mask}

process Local_Seeding_Mask {
    cpus 1

    input:
    set sid, file(wm), file(fa) from wm_fa_for_local_seeding_mask

    output:
    set sid, "${sid}__local_seeding_mask.nii.gz" into tracking_seeding_mask_for_local

    when:
        params.run_local_tracking

    script:
    if (params.local_seeding_mask_type == "wm")
        """
        mv $wm ${sid}__local_seeding_mask.nii.gz
        """
    else if (params.local_seeding_mask_type == "fa")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrcalc $fa $params.local_fa_seeding_mask_threshold -ge ${sid}__local_seeding_mask.nii.gz -datatype uint8
        """
}

fodf_for_local_tracking
    .join(tracking_mask_for_local)
    .join(tracking_seeding_mask_for_local)
    .set{fodf_maps_for_local_tracking}

process Local_Tracking {
    cpus 2

    input:
    set sid, file(fodf), file(tracking_mask), file(seed)\
        from fodf_maps_for_local_tracking
    each curr_seed from local_random_seed

    output:
    file "${sid}__local_tracking_${params.local_algo}_${params.local_seeding_mask_type}_seeding_${params.local_tracking_mask_type}_mask_seed_${curr_seed}.trk"

    when:
        params.run_local_tracking

    script:
    compress =\
        params.local_compress_streamlines ? '--compress ' + params.local_compress_value : ''
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_local_tracking.py $fodf $seed $tracking_mask\
            ${sid}__local_tracking_${params.local_algo}_${params.local_seeding_mask_type}_seeding_${params.local_tracking_mask_type}_mask_seed_${curr_seed}.trk\
            --algo $params.local_algo --$params.local_seeding $params.local_nbr_seeds\
            --seed $curr_seed --step $params.local_step --theta $params.local_theta\
            --sfthres $params.local_sfthres --min_length $params.local_min_len\
            --max_length $params.local_max_len $compress --sh_basis $params.basis
        """
}
