function [ll_list,ll_list_descreet,ll_list_rate] = general_loglikelihood(dist_obj,params,data,truncation_value)
%GENERAL_LOGLIKELIHOOD Summary of this function goes here
%   Detailed explanation goes here
ll_list = zeros(1,length(data));
ll_list_descreet = zeros(1,length(data));
ll_list_rate = zeros(1,length(data));
samples_per_cell = funoncellarray1input(data,@length);
samples_per_cell(isnan(samples_per_cell))=0;
for i =1:length(data)
    if isempty(data{i})
        per_rate(i)=0;
        ll_list(i) = nan;
        ll_list_descreet(i) = nan;
        continue;
    end
    
    data_scaled =data{i};
    data_upper = dbm2linear((linear2dbm(data_scaled)+.5));
    data_lower = dbm2linear((linear2dbm(data_scaled)-.5));
    trunc_val = truncation_value(i);
    trunc_val = dbm2linear((linear2dbm(trunc_val)-.5));
    pdf_nak_estimated = dist_obj.dist_handle(params(i,:));
    per_rate  = cdf(pdf_nak_estimated,trunc_val);
    model_pdf_samples = pdf(pdf_nak_estimated,data_scaled);
    model_pdf_samples = model_pdf_samples*(1/(1-per_rate));
    model_pdf_samples(model_pdf_samples<eps)=eps;
    model_pdf_samples(data_scaled<trunc_val)=eps;
    
    
    model_pmf_samples = cdf(pdf_nak_estimated,data_upper)-cdf(pdf_nak_estimated,data_lower);
    model_pmf_samples = model_pmf_samples*(1/(1-per_rate));
    model_pmf_samples(model_pmf_samples<eps)=eps;
    model_pmf_samples(data_scaled<trunc_val)=eps;
    
    ll_list(i) = -mean(log2(model_pdf_samples));
    ll_list_descreet(i) = -mean(log2(model_pmf_samples));
    
    
    data_scaled =pdf_nak_estimated.random(1e4,1);
    data_upper = dbm2linear((linear2dbm(data_scaled)+.5));
    data_lower = dbm2linear((linear2dbm(data_scaled)-.5));
    model_pmf_samples = cdf(pdf_nak_estimated,data_upper)-cdf(pdf_nak_estimated,data_lower);
    model_pmf_samples = model_pmf_samples*(1/(1-per_rate));
    model_pmf_samples(model_pmf_samples<eps)=eps;
    model_pmf_samples =model_pmf_samples(data_lower>=trunc_val);
    
    ll_list_desc_ehsan_entropy = -mean(log2(model_pmf_samples));
    ll_list_rate(i) = (ll_list_descreet(i)- ll_list_desc_ehsan_entropy)./ll_list_desc_ehsan_entropy;
    
    
end

end

