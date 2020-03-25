# covidm_interventions.R
# setting up intervention scenarios

# set up interventions
# TODO remove reliance on cm_translate_parameters, for efficiency?
cm_iv_build = function(parameters)
{
    p = cm_translate_parameters(parameters)
    data.table(date = (ymd(p$date0) + p$time0):(ymd(p$date0) + p$time1),
        contact = list(c(1, 1, 1, 1)));
}

cm_iv_checkset = function(iv, what, x = 1)
{
    if (!(what %in% names(iv))) {
        iv[, (what) := x];
    }
}

# set up school breaks
cm_iv_school_breaks = function(iv, ymd_break_start, ymd_break_end, sf = 0)
{
    cm_iv_contact(iv, ymd_break_start, ymd_break_end, c(1, 1, sf, 1));
}

# set up other interventions
# TODO add a trace component
# TODO min option instead of *
cm_iv_contact = function(iv, ymd_iv_first_day, ymd_iv_last_day, cf)
{
    # Check column is present
    cm_iv_checkset(iv, "contact");
    
    # Check start & end dates are same length
    if (length(ymd_iv_first_day) != length(ymd_iv_last_day)) {
        stop("length(ymd_iv_first_day) != length(ymd_iv_last_day)");
    }
    
    for (i in 1:length(ymd_iv_first_day)) {
        
        # Convert start and end dates
        t0 = ymd(ymd_iv_first_day[i]);
        t1 = ymd(ymd_iv_last_day[i]);
        
        # Apply intervention
        iv[date >= t0 & date <= t1, contact := lapply(contact, function(x) x * cf)];
    }
}

# apply interventions to a parameter set
# TODO merge with existing schedule
cm_iv_apply = function(parameters, iv)
{
    # Make schedule of interventions
    schedule = list();
    for (ivr in seq_len(nrow(iv)))
    {
        # Skip this row if nothing has changed
        if (ivr == 1) {
        } else if (any(unlist(iv[ivr, .SD, .SDcols = -"date"]) != unlist(iv[ivr - 1, .SD, .SDcols = -"date"]))) {
        } else next;
        
        # Set up changes
        change = list(t = as.numeric(as_date(iv[ivr, date]) - ymd(parameters$date0)));
        for (iv_col in 2:ncol(iv))
        {
            change[[length(change) + 1]] = iv[ivr, ..iv_col][[1]][[1]];
            names(change)[iv_col] = names(iv)[iv_col];
        }
        
        # Add to schedule
        schedule[[length(schedule) + 1]] = change;
    }

    # Apply to each population
    for (pi in 1:length(parameters$pop))
    {
        if (length(parameters$pop[[pi]]$schedule) == 0) {
            parameters$pop[[pi]]$schedule = schedule;
        } else {
            stop(paste0("Schedule already set for population ", pi, ". (Merging of schedules not yet implemented.)"));
        }
    }
    
    return (parameters);
}
