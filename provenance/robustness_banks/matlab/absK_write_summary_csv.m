function absK_write_summary_csv(summary,outcsv)
T = struct2table(summary,'AsArray',true);
writetable(T,outcsv);
end
