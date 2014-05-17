while true
	stmt = input('py> ', 's');
	if strcmp(stmt, 'exit')
		break
	else
		try
			py('eval', stmt);
		catch
		end
	end
end