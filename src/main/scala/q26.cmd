# commands to clean tables

cat item.dat | cut -d, -f1,10,13|  sed 's/[.!@#$%^&*()-]//g' | sed -e 's,_,,g' | tr abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ 012345678901234567890123456789  | cut -d ' '\
 -f1  > item_sanitized.csv
 
cat store_sales.dat | cut -d, -f3,4 > store_sales_sanitized.csv 
