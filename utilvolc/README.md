# Back Trajectories
  For the back trajectory analysis the steps are.

   *  create a csv file from the observations with the starting locations of the backwards trajectories.
      * This can be done with functions in utiltraj.py
         * trajectory_input_csv
      * Events class in volcat_event.py also has a method write_back_trajectory_csv
      * the csv file has columns 'lat','lon','mass','height','area' 
      * sometimes the csv file may have column 'heightI'. This is because the SO2 data for Hunga Tonga contains data that is interpolated. 
   * Run ash_main.py with 'backtrajectoryfromobs' as the runflag in the config file to run hyts_std.
      * generator functions in the trajectory_generators.py will read the csv file. 
        * the generator function is used in CollectTrajectory class.
        * currently the timegenerate_traj_from_obsdf generator is used. 
          It will creates an individual control file for each observation point in the csv file. The control file will run multiple heights from
          the same lat, lon point. Thus each HYSPLIT run represents one observation point and all the possible heights. The different heights
          are calculated in the generate_height_traj_from_series generator which takes a minimum height, maximum height, and height resolution as inputs.
        * In this configuration one HYSPLIT run per observation point is done.
        * Alternate configurations are to 
          * run all trajectories that start at the same time together.
          * do all the runs individually (e.g. only one trajectory per run)
   * After the runs are completed the combine_traj function in utiltraj.py may be used to read the information from all the tdump outputs into one dataframe.

   * this data frame may then be used in TrajAnalysis class to create various plots.

   * this dataframe may also be used by some functions in utiltraj.py to create emit-times files for data insertion.




# TCM


   utiltcm
       classes
       ParametersIn

       InverseDat - 
           InvEstimatedEmissions
           InverseOut2Dat
   

   volcinverse.py
       classes -
       InvDirectory
       InversionEns
       RunInversion
          uses VolcatHysplit paired data class.
      

   volctcm

       FUNCTIONS - 
           remove_near_clear_sky(avg, window)

       TCM class.
           represents the TCM.
           attributes:
               _columns, _tcm
               n_ctrl
               tcm_name
               output
               tag
               latlist
               longlist
               tcm_lat
               tcm_lon
           properties - tcm (readonly), columns
           methods
               plot
               run
               make_tcm
               make_tcm_mult
               make_tcm_names
               write 
           

   volcpaired
       VolcatHysplit class
 

# MODEL+VOLCAT EVALUATION

   volcpaired.py
       VolcatHysplit class


# VOLCAT

    * volcat.py for reading VOLCAT netcdf event files.
       * DONE - TODO - _get_time function needs checking / updating
       * TODO - find_volcat function may need simplifying / splitting.
       * TODO - might make some of the functions into a class with an interface.


    * volcat_files.py 
       * TODO - may be some unused functions in here.

    * volcat_plots.py
       * TODO - the box plot and cdf plotting functionality needs updating.
              

    * volcat_2dplots.py


    * TODO - need more functions/methods for identifying/combining volcat files that are at same time
           but have different id numbers. 
           this has been started in volcat_events with the OneEventTime class.

   * volcat_event.py

        FUNCTIONS 
            create_event_from_fnames
            create_event 

        Events class

          Attributes:
              self.eventdf - EventDF class

          Properties
              edf - dataframe of self.eventdf 

          Probably has too many methods!

        * EventDisplay class (helper class)
        * EventStatus  class (helper class)
        * EventDF      class (helper class)
        * OneEventTime class (helper class)


# Model data
   * model_event.py contains classes for helping with QVA creation.

        * ModelEvent class 
             * DataInsertion outputs
             *  Model runs for input into inversion algorithm.
             * Model runs which are result of emissions with inversion algorithm.
             * Model runs which are the result of 'traditional' run.

        * ModelForecast class
             Instances of the forecast class are created by the ModelEvent class.
             the forecast class should also be able to ingest the VOLCAT data.

   ensemble_polygons contains classes for taking model data and creating polygons

 
   * qvainterface.py
   TODO - just started writing this.




# INTERFACES

    ashapp/ashruninterface  classes for running HYSPLIT model, postprocessing output and creating graphics.

           DONE - the dispersion model runs will run automatically.
           TODO - possibly add some type hinting.
    
    utilvolc/inversioninterface 


# DATA INSERTION
    
    ashapp/DImain.py 
    ashapp/make_volcat_database.py  run this to create a csv file from the sumdf file created by get_summary_file_df.





