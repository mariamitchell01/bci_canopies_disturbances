@ECHO OFF
ECHO.
ECHO Your computer name is:
hostname
ECHO.
echo.
ipconfig
ECHO.
getmac
ECHO.
echo.
systeminfo | findstr /B /C:"OS Name" /C:"OS Version"
echo.
echo Your computer name name is:
hostname
echo.
echo Service Desk
echo Phone: x34000 or 202-633-4000
echo Email: OCIOservicedesk@si.edu
echo https://smithsonianprod.servicenowservices.com/si
echo.
PAUSE