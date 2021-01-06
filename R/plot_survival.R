#' Plot a Kaplan-Meier Curve
#' @param dbms one of the following supported Database Management Systems: c("oracle","postgresql","redshift","sql server","pdw", "netezza","bigquery","sqlite")
#' @param user username for OMOP v5.0 instance
#' @param password password for OMOP v5.0 instance
#' @param server OMOP v5.0 server
#' @param port Port.

plot_survival <-
        function(dbms = c("oracle","postgresql","redshift","sql server","pdw", "netezza","bigquery","sqlite"),
                 user = NULL,
                 password = NULL,
                 server,
                 port,
                 omop_cdm_schema = "omop_cdm") {

        con_details <-
        DatabaseConnector::createConnectionDetails(dbms = dbms,
                                                   user = user,
                                                   password = password,
                                                   server = server,
                                                   port = port)

        con <- DatabaseConnector::connect(connectionDetails = con_details)

        rendered_sql <-
                SqlRender::loadRenderTranslateSql(sqlFilename = "time_dx_to_survival.sql",
                                                  packageName = "oncoPlot",
                                                  dbms = dbms,
                                                  cdmSchema = schema)

        rendered_sql <-
                SqlRender::render(
                        sql =
                                "
                                -- Survival from diagnosis, record level
                                select
                                ed.person_id,
                                c.concept_name as cohort_cols,
                                meas.is_metastatic,
                                case WHEN dth.death_datetime is not null then 1 else 0 end as event_col,
                                --DATEDIFF(month, ed.episode_start_datetime, coalesce(dth.death_datetime, op.last_followup_date)) as survival_time_col
                                ed.episode_start_datetime as start_date,
                                coalesce(dth.death_datetime, op.last_followup_date) as end_date
                                from @cdmSchema.episode ed
                                left join @cdmSchema.person p
                                on ed.episode_concept_id = 32528 -- first disease occurrence
                                and ed.person_id = p.person_id
                                left join
                                ( select distinct person_id , MAX(death_datetime) death_datetime
                                from @cdmSchema.death
                                group by person_id
                                )  dth
                                on p.person_id = dth.person_id
                                left join
                                (select max(observation_period_end_date) as last_followup_date, person_id
                                from @cdmSchema.observation_period
                                group by person_id
                                )  op
                                on ed.person_id = op.person_id
                                join @cdmSchema.concept_relationship cr1
                                on ed.episode_object_concept_id = cr1.concept_id_2
                                and cr1.relationship_id = 'Maps to'
                                left join @cdmSchema.concept_relationship cr2
                                on cr1.concept_id_1 = cr2.concept_id_1
                                and cr2.relationship_id = 'ICDO to Schema'
                                join @cdmSchema.concept c
                                on cr2.concept_id_2 = c.concept_id
                                left join
                                (
                                  SELECT modifier_of_event_id
                                		,CASE
                                			WHEN metastatic > 0 THEN 'Metastatic'
                                			WHEN metastatic = 0 AND nonmeta > 0 THEN 'Non-Metastatic'
                                			ELSE 'Unknown'
                                		END as is_metastatic
                                 FROM
                                 (
                                  SELECT modifier_of_event_id
                                		,COUNT(CASE WHEN value_source_value IN ('1','2') THEN 1 END) metastatic
                                		,COUNT(CASE WHEN value_source_value IN('0')  THEN 1 END ) nonmeta
                                 FROM @cdmSchema.measurement
                                 WHERE modifier_of_field_concept_id = 1000000003 -- 'epsiode.episode_id'
                                 AND measurement_concept_id
                                 IN (

                                 35918581 -- Mets at DX-Bone
                                ,35918692 -- Mets at DX-Brain
                                ,35918491 -- Mets at DX-Distant LN
                                ,35918290 -- Mets at DX-Liver
                                ,35918559 -- Mets at DX-Lung
                                ,35918527 -- Mets at DX-Other
                                 )
                                 GROUP BY modifier_of_event_id
                                 ) x
                                 ) meas
                                 on ed.episode_id = meas.modifier_of_event_id
                        ",
                        cdmSchema = omop_cdm_schema
                )

        rendered_sql <-
        SqlRender::translate(sql = rendered_sql,
                             targetDialect = dbms)

        dataframe <- DatabaseConnector::dbGetQuery(con, statement = rendered_sql)

        DatabaseConnector::dbDisconnect(conn = con)

        # trim string names
        dataframe$cohort_cols <- sub("[(].*", "", dataframe$cohort_cols)

        # filter to only the top 10 most prevalent
        top10 <- names(sort(table(dataframe$cohort_cols), decreasing = TRUE)[1:10])
        dataframe <- dataframe[dataframe$cohort_cols %in% top10,]

        dataframe$event_col <- as.numeric(dataframe$event_col)

        dataframe$cohort_cols <- as.factor(dataframe$cohort_cols)

        dataframe$timediff <- difftime(as.Date(dataframe$end_date),as.Date(dataframe$start_date), units = "days")/30

        survival_object <<- try_catch_error_as_na(survival::Surv(time = dataframe$timediff,
                                                       event = dataframe$event_col,
                                                       type = "right"))

        if (is.vector(survival_object)) {
                cat("\n\tError: survival_time and/or event_occurred not in correct format. Please check and try again.\n")
        } else {
                km_fit_01 <- try_catch_error_as_na(survival::survfit(survival_object ~ cohort_cols,
                                                           data = dataframe))
                if ((length(km_fit_01) == 1) & any(is.na(km_fit_01))) {
                        cat("\n\tError: cohort_object and/or dataframe not in correct format. Please check and try again.\n")

                } else {
                        medsurv <- survminer::surv_median(km_fit_01)

                        OUTPUT <- survminer::ggsurvplot(km_fit_01,
                                             data = dataframe,
                                             xscale = 12,
                                             break.x.by = 6,
                                             legend = c(0.8, 0.9),
                                             surv.median.line = "hv",
                                             legend.title = "Cohort",
                                             legend.labs = levels(dataframe %>%
                                                                          select(cohort_cols) %>% unlist())) + ggplot2::xlab("Survival Time (Years)") + ggplot2::ylab("Survival Probability") + ggplot2::ggtitle("Kaplan-Meier Curves")

                        OUTPUT$plot + ggplot2::annotate("text",
                                                        x = medsurv$median + 2,
                                                        y = (1:nrow(medsurv))/20,
                                                        label = round(medsurv$median/12, 2),
                                                        parse = TRUE)
                }

        }
}


