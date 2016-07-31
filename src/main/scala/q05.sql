SELECT
  --wcs_user_sk,
  clicks_in_category,
  CASE WHEN cd_education_status IN (-education) THEN 1 ELSE 0 END AS college_education,
  CASE WHEN cd_gender = (-gender) THEN 1 ELSE 0 END AS male,
  clicks_in_1,
  clicks_in_2,
  clicks_in_3,
  clicks_in_4,
  clicks_in_5,
  clicks_in_6,
  clicks_in_7
FROM( 
  SELECT 
    wcs_user_sk,
    SUM( CASE WHEN i_category = (-category) THEN 1 ELSE 0 END) AS clicks_in_category,
    SUM( CASE WHEN i_category_id = 1 THEN 1 ELSE 0 END) AS clicks_in_1,
    SUM( CASE WHEN i_category_id = 2 THEN 1 ELSE 0 END) AS clicks_in_2,
    SUM( CASE WHEN i_category_id = 3 THEN 1 ELSE 0 END) AS clicks_in_3,
    SUM( CASE WHEN i_category_id = 4 THEN 1 ELSE 0 END) AS clicks_in_4,
    SUM( CASE WHEN i_category_id = 5 THEN 1 ELSE 0 END) AS clicks_in_5,
    SUM( CASE WHEN i_category_id = 6 THEN 1 ELSE 0 END) AS clicks_in_6,
    SUM( CASE WHEN i_category_id = 7 THEN 1 ELSE 0 END) AS clicks_in_7
  FROM web_clickstreams
  INNER JOIN item it ON (wcs_item_sk = i_item_sk
                     AND wcs_user_sk IS NOT NULL)
  GROUP BY  wcs_user_sk
)q05_user_clicks_in_cat
INNER JOIN customer ct ON wcs_user_sk = c_customer_sk
INNER JOIN customer_demographics ON c_current_cdemo_sk = cd_demo_sk
;
