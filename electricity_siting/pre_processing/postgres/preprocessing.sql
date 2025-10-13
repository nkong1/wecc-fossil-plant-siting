-- Coal without CCS
ALTER TABLE public."coal_igcc — wecc_gridcerf_coal_igcc"
  ADD COLUMN IF NOT EXISTS dist_to_surface_flow_meters FLOAT;

CREATE TEMP TABLE wecc_coal_water_dists AS
SELECT
  g.fid,
  MIN(ST_Distance(g.geom, p.geom)) AS dist_to_surface_flow_meters
FROM public."coal_igcc — wecc_gridcerf_coal_igcc" AS g
JOIN public.clipped_greaterthan95mgd AS p
  ON ST_DWithin(g.geom, p.geom, 30000)
GROUP BY g.fid;

UPDATE public."coal_igcc — wecc_gridcerf_coal_igcc" AS main
SET dist_to_surface_flow_meters = d.dist_to_surface_flow_meters
FROM wecc_coal_water_dists AS d
WHERE main.fid = d.fid;


-- Coal with CCS
ALTER TABLE public."coal_igcc_ccs — wecc_gridcerf_coal_igcc_ccs"
  ADD COLUMN IF NOT EXISTS dist_to_surface_flow_meters FLOAT;

CREATE TEMP TABLE wecc_coal_ccs_water_dists3 AS
SELECT
  g.fid,
  MIN(ST_Distance(g.geom, p.geom)) AS dist_to_surface_flow_meters
FROM public."coal_igcc_ccs — wecc_gridcerf_coal_igcc_ccs" AS g
JOIN public.clipped_greaterthan145mgd AS p
  ON ST_DWithin(g.geom, p.geom, 30000)
GROUP BY g.fid;

UPDATE public."coal_igcc_ccs — wecc_gridcerf_coal_igcc_ccs" AS main
SET dist_to_surface_flow_meters = d.dist_to_surface_flow_meters
FROM wecc_coal_ccs_water_dists AS d
WHERE main.fid = d.fid;


-- Gas without CCS
ALTER TABLE public."gas_cc — wecc_gridcerf_gas_cc"
  ADD COLUMN IF NOT EXISTS dist_to_surface_flow_meters FLOAT;

CREATE TEMP TABLE wecc_gas_water_dists AS
SELECT
  g.fid,
  MIN(ST_Distance(g.geom, p.geom)) AS dist_to_surface_flow_meters
FROM public."gas_cc — wecc_gridcerf_gas_cc" AS g
JOIN public.clipped_greaterthan70mgd AS p
  ON ST_DWithin(g.geom, p.geom, 30000)
GROUP BY g.fid;

UPDATE public."gas_cc — wecc_gridcerf_gas_cc" AS main
SET dist_to_surface_flow_meters = d.dist_to_surface_flow_meters
FROM wecc_gas_water_dists AS d
WHERE main.fid = d.fid;


-- Gas with CCS
ALTER TABLE public."gas_cc_ccs — wecc_gridcerf_gas_cc_ccs"
  ADD COLUMN IF NOT EXISTS dist_to_surface_flow_meters FLOAT;

CREATE TEMP TABLE wecc_gas_ccs_water_dists AS
SELECT
  g.fid,
  MIN(ST_Distance(g.geom, p.geom)) AS dist_to_surface_flow_meters
FROM public."gas_cc_ccs — wecc_gridcerf_gas_cc_ccs" AS g
JOIN public.clipped_greaterthan120mgd AS p
  ON ST_DWithin(g.geom, p.geom, 30000)
GROUP BY g.fid;

UPDATE public."gas_cc_ccs — wecc_gridcerf_gas_cc_ccs" AS main
SET dist_to_surface_flow_meters = d.dist_to_surface_flow_meters
FROM wecc_gas_ccs_water_dists AS d
WHERE main.fid = d.fid; 

