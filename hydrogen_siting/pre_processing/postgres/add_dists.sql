
------------------------------------------------------------------------------
-- Assign all layers lat/lon coords:
------------------------------------------------------------------------------

--biogas smr with ccs
ALTER TABLE public.bio_smr_ccs
ADD COLUMN latitude double precision, 
ADD COLUMN longitude double precision;

UPDATE public.bio_smr_ccs
SET
  longitude = ST_X(ST_CENTROID(ST_TRANSFORM(geom, 4326))),
  latitude = ST_Y(ST_CENTROID(ST_TRANSFORM(geom, 4326)));


--biogas smr
ALTER TABLE public.bio_smr
ADD COLUMN latitude double precision, 
ADD COLUMN longitude double precision;

UPDATE public.bio_smr
SET
  longitude = ST_X(ST_CENTROID(ST_TRANSFORM(geom, 4326))),
  latitude = ST_Y(ST_CENTROID(ST_TRANSFORM(geom, 4326)));


--biomass gasification
ALTER TABLE public.biomass_gas
ADD COLUMN latitude double precision, 
ADD COLUMN longitude double precision;

UPDATE public.biomass_gas
SET
  longitude = ST_X(ST_CENTROID(ST_TRANSFORM(geom, 4326))),
  latitude = ST_Y(ST_CENTROID(ST_TRANSFORM(geom, 4326)));


--coal gasification with ccs
ALTER TABLE public.coal_gas_ccs
ADD COLUMN latitude double precision, 
ADD COLUMN longitude double precision;

UPDATE public.coal_gas_ccs
SET
  longitude = ST_X(ST_CENTROID(ST_TRANSFORM(geom, 4326))),
  latitude = ST_Y(ST_CENTROID(ST_TRANSFORM(geom, 4326)));


--coal gasification 
ALTER TABLE public.coal_gas
ADD COLUMN latitude double precision, 
ADD COLUMN longitude double precision;

UPDATE public.coal_gas
SET
  longitude = ST_X(ST_CENTROID(ST_TRANSFORM(geom, 4326))),
  latitude = ST_Y(ST_CENTROID(ST_TRANSFORM(geom, 4326)));


--natural gas smr
ALTER TABLE public.gas_smr
ADD COLUMN latitude double precision, 
ADD COLUMN longitude double precision;

UPDATE public.gas_smr
SET
  longitude = ST_X(ST_CENTROID(ST_TRANSFORM(geom, 4326))),
  latitude = ST_Y(ST_CENTROID(ST_TRANSFORM(geom, 4326)));


------------------------------------------------------------------------------
-- Create and fill in columns for closest distance to feedstock source:
------------------------------------------------------------------------------

--biogas smr with ccs
CREATE TEMP TABLE biogas_ccs_dists_new AS
SELECT
  g.fid,
  ST_Distance(g.geom, p.geom) AS new_dist
FROM public.bio_smr_ccs AS g
JOIN LATERAL (
	SELECT geom
	  FROM (
	    SELECT geom FROM public."eia_biogas_plants — eia_biogas_reprojected"
	    UNION ALL 
	    SELECT geom FROM public."livestock-digesters"
		UNION ALL 
	    SELECT geom FROM public.landfills
		UNION ALL 
	    SELECT geom FROM public.wastewater_plants
	) AS p
	  WHERE ST_DWithin(g.geom, p.geom, 10000)
  ORDER BY g.geom <-> p.geom
  LIMIT 1
) AS p ON true;


ALTER TABLE public.bio_smr_ccs
ADD COLUMN dist_to_feedstock_meters double precision;


UPDATE public.bio_smr_ccs AS main
SET dist_to_feedstock_meters = new.new_dist
FROM biogas_ccs_dists_new AS new
WHERE main.fid = new.fid;


--biogas smr 
CREATE TEMP TABLE biogas_dists_new AS
SELECT
  g.fid,
  ST_Distance(g.geom, p.geom) AS new_dist
FROM public.bio_smr AS g
JOIN LATERAL (
	SELECT geom
	  FROM (
	    SELECT geom FROM public."eia_biogas_plants — eia_biogas_reprojected"
	    UNION ALL 
	    SELECT geom FROM public."livestock-digesters"
		UNION ALL 
	    SELECT geom FROM public.landfills
		UNION ALL 
	    SELECT geom FROM public.wastewater_plants
	) AS p
	  WHERE ST_DWithin(g.geom, p.geom, 10000)
  ORDER BY g.geom <-> p.geom
  LIMIT 1
) AS p ON true;


ALTER TABLE public.bio_smr
ADD COLUMN dist_to_feedstock_meters double precision;

UPDATE public.bio_smr AS main
SET dist_to_feedstock_meters = new.new_dist
FROM biogas_dists_new AS new
WHERE main.fid = new.fid;


--biomass gasification
CREATE TEMP TABLE biomass_dists AS
SELECT
  g.fid,
  ST_Distance(g.geom, p.geom) AS new_dist
FROM public.biomass_gas AS g
JOIN LATERAL (
	SELECT geom 
	  FROM (
	    SELECT geom FROM public."eia_biomass — biomass"
	) AS p
	  WHERE ST_DWithin(g.geom, p.geom, 5000)
  ORDER BY g.geom <-> p.geom
  LIMIT 1
) AS p ON true;

ALTER TABLE public.biomass_gas
ADD COLUMN dist_to_feedstock_meters double precision;

UPDATE public.biomass_gas AS main
SET dist_to_feedstock_meters = new.new_dist
FROM biomass_dists AS new
WHERE main.fid = new.fid;


--coal gasification with ccs
CREATE TEMP TABLE coal_ccs_dists AS
SELECT
  g.fid,
  ST_Distance(g.geom, p.geom) AS new_dist
FROM public.coal_gas_ccs AS g
JOIN LATERAL (
	SELECT geom
	  FROM (
	    SELECT geom FROM public.coal_mines
		UNION ALL 
		SELECT geom FROM public."coal_plants — eia_860_2024_er_gen_w_coords20250730t183445z100"
	) AS p
  ORDER BY g.geom <-> p.geom
  LIMIT 1
) AS p ON true;

ALTER TABLE public.coal_gas_ccs
ADD COLUMN dist_to_feedstock_meters double precision;

UPDATE public.coal_gas_ccs AS main
SET dist_to_feedstock_meters = new.new_dist
FROM coal_ccs_dists AS new
WHERE main.fid = new.fid;



--coal gasification
CREATE TEMP TABLE coal_dists AS
SELECT
  g.fid,
  ST_Distance(g.geom, p.geom) AS new_dist
FROM public.coal_gas AS g
JOIN LATERAL (
	SELECT geom
	  FROM (
	    SELECT geom FROM public.coal_mines
		UNION ALL 
		SELECT geom FROM public."coal_plants — eia_860_2024_er_gen_w_coords20250730t183445z100"
	) AS p
  ORDER BY g.geom <-> p.geom
  LIMIT 1
) AS p ON true;

ALTER TABLE public.coal_gas
ADD COLUMN dist_to_feedstock_meters double precision;

UPDATE public.coal_gas AS main
SET dist_to_feedstock_meters = new.new_dist
FROM coal_dists AS new
WHERE main.fid = new.fid;


--natural gas smr
CREATE TEMP TABLE gas_smr_dists AS
SELECT
  g.fid,
  ST_Distance(g.geom, p.geom) AS new_dist
FROM public.gas_smr AS g
JOIN LATERAL (
	SELECT geom
	  FROM (
	    SELECT geom FROM public.clipped_pipelines
	) AS p
  ORDER BY g.geom <-> p.geom
  LIMIT 1
) AS p ON true;

ALTER TABLE public.gas_smr
ADD COLUMN dist_to_feedstock_meters double precision;

UPDATE public.gas_smr AS main
SET dist_to_feedstock_meters = new.new_dist
FROM gas_smr_dists AS new
WHERE main.fid = new.fid;


------------------------------------------------------------------------------
-- Create and fill in columns for closest distance to electric substation:
------------------------------------------------------------------------------
--biogas smr with ccs
CREATE TEMP TABLE biogas_ccs_substation_dists AS
SELECT
  g.fid,
  ST_Distance(g.geom, p.geom) AS new_dist
FROM public.bio_smr_ccs AS g
JOIN LATERAL (
	SELECT geom
	  FROM (
	    SELECT geom FROM public.substations_clipped
	) AS p
  ORDER BY g.geom <-> p.geom
  LIMIT 1
) AS p ON true;

ALTER TABLE public.bio_smr_ccs
ADD COLUMN dist_to_substation_meters double precision;

UPDATE public.bio_smr_ccs AS main
SET dist_to_substation_meters = new.new_dist
FROM biogas_ccs_substation_dists AS new
WHERE main.fid = new.fid;


--biogas smr 
CREATE TEMP TABLE biogas_substation_dists AS
SELECT
  g.fid,
  ST_Distance(g.geom, p.geom) AS new_dist
FROM public.bio_smr AS g
JOIN LATERAL (
	SELECT geom
	  FROM (
	    SELECT geom FROM public.substations_clipped
	) AS p
  ORDER BY g.geom <-> p.geom
  LIMIT 1
) AS p ON true;

ALTER TABLE public.bio_smr
ADD COLUMN dist_to_substation_meters double precision;

UPDATE public.bio_smr AS main
SET dist_to_substation_meters = new.new_dist
FROM biogas_substation_dists AS new
WHERE main.fid = new.fid;


--biomass gasification
CREATE TEMP TABLE biomass_substation_dists AS
SELECT
  g.fid,
  ST_Distance(g.geom, p.geom) AS new_dist
FROM public.biomass_gas AS g
JOIN LATERAL (
	SELECT geom 
	  FROM (
	    SELECT geom FROM public.substations_clipped
	) AS p
  ORDER BY g.geom <-> p.geom
  LIMIT 1
) AS p ON true;

ALTER TABLE public.biomass_gas
ADD COLUMN dist_to_substation_meters double precision;

UPDATE public.biomass_gas AS main
SET dist_to_substation_meters = new.new_dist
FROM biomass_substation_dists AS new
WHERE main.fid = new.fid;


--coal gasification with ccs
CREATE TEMP TABLE coal_ccs_substation_dists AS
SELECT
  g.fid,
  ST_Distance(g.geom, p.geom) AS new_dist
FROM public.coal_gas_ccs AS g
JOIN LATERAL (
	SELECT geom
	  FROM (
	    SELECT geom FROM public.substations_clipped
	) AS p
  ORDER BY g.geom <-> p.geom
  LIMIT 1
) AS p ON true;

ALTER TABLE public.coal_gas_ccs
ADD COLUMN dist_to_substation_meters double precision;

UPDATE public.coal_gas_ccs AS main
SET dist_to_substation_meters = new.new_dist
FROM coal_ccs_substation_dists AS new
WHERE main.fid = new.fid;



--coal gasification
CREATE TEMP TABLE coal_substation_dists AS
SELECT
  g.fid,
  ST_Distance(g.geom, p.geom) AS new_dist
FROM public.coal_gas AS g
JOIN LATERAL (
	SELECT geom
	  FROM (
	    SELECT geom FROM public.substations_clipped
	) AS p
  ORDER BY g.geom <-> p.geom
  LIMIT 1
) AS p ON true;

ALTER TABLE public.coal_gas
ADD COLUMN dist_to_substation_meters double precision;

UPDATE public.coal_gas AS main
SET dist_to_substation_meters = new.new_dist
FROM coal_substation_dists AS new
WHERE main.fid = new.fid;


--natural gas smr
CREATE TEMP TABLE gas_smr_substation_dists AS
SELECT
  g.fid,
  ST_Distance(g.geom, p.geom) AS new_dist
FROM public.gas_smr AS g
JOIN LATERAL (
	SELECT geom
	  FROM (
	    SELECT geom FROM public.substations_clipped
	) AS p
  ORDER BY g.geom <-> p.geom
  LIMIT 1
) AS p ON true;

ALTER TABLE public.gas_smr
ADD COLUMN dist_to_substation_meters double precision;

UPDATE public.gas_smr AS main
SET dist_to_substation_meters = new.new_dist
FROM gas_smr_substation_dists AS new
WHERE main.fid = new.fid;
